#!/usr/bin/env python3
"""Build a quantification spec TSV from a YAML definition.

All biology lives in the YAML; this script only expands archetype rules into
the flat TSV that R reads at runtime. Myeloid coarse classes are always derived
from fine (Classify_Myeloid collapse in RIRA), not taken from YAML as source of
truth.

  python3 inst/scripts/build_quantify_spec.py quantify_human_immune.yml \\
    -o inst/extdata/quantify_human_immune_spec.tsv
  python3 inst/scripts/build_quantify_spec.py quantify_rhesus.yml \\
    -o inst/extdata/quantify_rhesus_spec.tsv \\
    --models-dir /path/to/RIRA/inst/models
"""
from __future__ import annotations

import argparse
import pickle
import re
from pathlib import Path

import yaml

COLUMNS = [
    "GroupingVariable", "TargetField", "SourceField", "MatchValue",
    "QuantificationType", "QuantificationSourceField", "QuantificationScoreType",
    "ScopeField", "ScopeMatchValue", "EffectorDifferentiationScoreField",
    "EffectorDifferentiationCutpointLow", "EffectorDifferentiationCutpointHigh",
    "SubsetPhenotypeOutputFieldName",
]

# Post-pkl labels from RIRA Classify_* (CellTypist.R); appended when --models-dir.
PIPELINE_META = ["Unassigned", "Ambiguous", "Unknown"]
IMMUNE_DERIVED = ["Non-Immune", "Unknown"]
T_AND_NK_DERIVED = ["Other"]

# Myeloid fine → coarse collapse (RIRA Classify_Myeloid); not authored in YAML.
MONOCYTE_FINE = {"CD14+ Monocytes", "CD16+ Monocytes", "Inflammatory Monocytes"}
DC_FINE = {"DC", "Mature DC"}


def load_quantify_config(path: Path) -> dict:
    """Load a quantify YAML config; shared by spec and role-map builders."""
    return yaml.safe_load(path.read_text())


def derive_coarse(fine_classes: list[str]) -> list[str]:
    """Collapse monocyte/DC fine labels; append pipeline meta for coerce parity."""
    coarse: set[str] = set()
    for fine in fine_classes:
        if fine in MONOCYTE_FINE:
            coarse.add("Monocytes")
        elif fine in DC_FINE:
            coarse.add("DC")
        else:
            coarse.add(fine)
    for label in PIPELINE_META:
        coarse.add(label)
    return sorted(coarse, key=str.lower)


def _load_pkl_classes(pkl_path: Path) -> list[str]:
    with pkl_path.open("rb") as handle:
        return list(pickle.load(handle)["Model"].classes_)


def _append_unique(base: list[str], extra: list[str]) -> list[str]:
    return base + [c for c in extra if c not in base]


def refresh_classes_from_models(cfg: dict, models_dir: Path) -> None:
    """Replace field class lists from pkl; in-memory only (TSV output, not YAML)."""
    meta = cfg.get("pipeline_meta", PIPELINE_META)
    for name, field in cfg.get("fields", {}).items():
        pkl_name = field.get("pkl")
        if not pkl_name:
            continue
        pkl_classes = _load_pkl_classes(models_dir / pkl_name)
        if name == "immune":
            field["classes"] = _append_unique(pkl_classes, IMMUNE_DERIVED)
        elif name == "tnk":
            field["classes"] = _append_unique(
                _append_unique(pkl_classes, T_AND_NK_DERIVED), meta
            )
        elif name == "myeloid_fine":
            field["classes"] = _append_unique(pkl_classes, meta)
        elif name == "myeloid_coarse":
            # Coarse follows pkl fine collapse, not YAML fine list.
            field["classes"] = derive_coarse(pkl_classes)
        else:
            field["classes"] = _append_unique(pkl_classes, meta)


def apply_myeloid_coarse_derivation(cfg: dict) -> None:
    """Overwrite myeloid_coarse.classes from myeloid_fine; YAML coarse list is not authoritative."""
    fields = cfg.get("fields", {})
    fine = fields.get("myeloid_fine")
    coarse = fields.get("myeloid_coarse")
    if not fine or not coarse:
        return
    coarse["classes"] = derive_coarse(list(fine["classes"]))


def sanitize(label: str) -> str:
    s = label.replace("+", "plus")
    return re.sub(r"_+", "_", re.sub(r"[^\w]+", "_", s)).strip("_")


class Spec:
    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.grouping = cfg["grouping"]
        self.delim = cfg.get("delimiter", "__")
        self.meta = cfg.get("pipeline_meta", [])
        self.fields = cfg["fields"]
        self.roles = cfg.get("roles", {})
        self.rows: list[list[str]] = []

    def column(self, field):        return self.fields[field]["column"]
    def prefix(self, field):        return self.fields[field]["prefix"]
    def parentmap(self, field):     return self.fields[field].get("parent", {})
    def parent_of(self, field, c):  return self.parentmap(field).get(c)
    def parent_col(self, field):    return self.column(self.fields[field]["parent_field"])

    def classes(self, field):
        base = list(self.fields[field]["classes"])
        return base + [m for m in self.meta if m not in base]

    def role(self, name):           return self.roles[name]

    def select(self, field, sel):
        classes = self.classes(field)
        cset = set(classes)
        if sel in (None, "all"):
            return classes
        if "class" in sel:
            return [sel["class"]]
        if "class_in" in sel:
            allowed = set(self.role(sel["class_in"]))
            return [c for c in classes if c in allowed]
        if "parent" in sel:
            return [c for c in classes if self.parent_of(field, c) == sel["parent"]]
        if "parent_in" in sel:
            allowed = set(self.role(sel["parent_in"]))
            return [c for c in classes if self.parent_of(field, c) in allowed]
        if "role" in sel:
            names = sel["role"] if isinstance(sel["role"], list) else [sel["role"]]
            out = []
            for n in names:
                for c in self.role(n):
                    if c in cset and c not in out:
                        out.append(c)
            return out
        raise ValueError(f"bad selector: {sel}")

    def match_pred(self, field, c, pred):
        if "class" in pred:     return c == pred["class"]
        if "class_in" in pred:  return c in set(self.role(pred["class_in"]))
        if "parent" in pred:    return self.parent_of(field, c) == pred["parent"]
        if "parent_in" in pred: return self.parent_of(field, c) in set(self.role(pred["parent_in"]))
        return False

    def scope(self, field, c, scope):
        if not scope:                return "", ""
        if scope == "by_parent":
            p = self.parent_of(field, c)
            return (self.parent_col(field), p) if p else ("", "")
        if scope == "by_self":       return self.column(field), c
        if "by_parent_in" in scope:
            p = self.parent_of(field, c)
            return (self.parent_col(field), p) if p in set(scope["by_parent_in"]) else ("", "")
        return self.column(scope["field"]), scope["value"]

    def target(self, field, *parts):  return self.delim.join([self.prefix(field), *parts])
    def emit(self, **kw):
        row = {"GroupingVariable": self.grouping, **kw}
        self.rows.append([str(row.get(c, "")) for c in COLUMNS])


def a_sum(s, r):
    f = r["field"]
    for c in s.select(f, r.get("select", "all")):
        sf, sv = s.scope(f, c, r.get("scope"))
        s.emit(TargetField=s.target(f, sanitize(c)), SourceField=s.column(f),
               MatchValue=c, QuantificationType="sum", ScopeField=sf, ScopeMatchValue=sv)

def a_score(s, r):
    f = r["field"]
    when = r.get("modules_when")
    for c in s.select(f, r.get("select", "all")):
        if when:
            mods = next((w["modules"] for w in when if s.match_pred(f, c, w["when"])), None)
            if mods is None:
                continue
        else:
            mods = r["modules"]
        sf, sv = s.scope(f, c, r.get("scope"))
        for suffix, ucell in s.cfg["modules"][mods] if isinstance(mods, str) else mods:
            s.emit(TargetField=s.target(f, sanitize(c), suffix), SourceField=s.column(f),
                   MatchValue=c, QuantificationType="score", QuantificationSourceField=ucell,
                   QuantificationScoreType="quantiles", ScopeField=sf, ScopeMatchValue=sv)

def a_eds(s, r):
    f, eds = r["field"], s.cfg["eds"]
    for c in s.select(f, r["select"]):
        sf, sv = s.scope(f, c, r.get("scope"))
        for ph in eds["phenotypes"]:
            s.emit(TargetField=s.target(f, sanitize(c), ph), SourceField=s.column(f),
                   MatchValue=c, QuantificationType="sum", ScopeField=sf, ScopeMatchValue=sv,
                   EffectorDifferentiationScoreField=eds["metadata_field"],
                   EffectorDifferentiationCutpointLow=eds["cutpoint_low"],
                   EffectorDifferentiationCutpointHigh=eds["cutpoint_high"],
                   SubsetPhenotypeOutputFieldName=ph)
        s.emit(TargetField=s.target(f, sanitize(c), "EDSScore"), SourceField=s.column(f),
               MatchValue=c, QuantificationType="score",
               QuantificationSourceField=eds["metadata_field"], QuantificationScoreType="quantiles",
               ScopeField=sf, ScopeMatchValue=sv)

def a_pct_positive(s, r):
    f, pc = r["field"], s.cfg["pct_positive"]
    order = s.role(pc["order_by"])
    cset = set(s.classes(f))
    for gene, spec in pc["genes"].items():
        names = pc["default_roles"] if spec == "default" else spec
        wanted = {c for n in names for c in s.role(n)}
        for c in order:
            if c in wanted and c in cset:
                sf, sv = s.scope(f, c, r.get("scope"))
                s.emit(TargetField=s.target(f, sanitize(c), "PctPositive", gene),
                       SourceField=s.column(f), MatchValue=c, QuantificationType="pct_positive",
                       QuantificationSourceField=gene, ScopeField=sf, ScopeMatchValue=sv)

def a_diversity(s, r):
    f = r["field"]
    for c in s.select(f, r["select"]):
        sf, sv = s.scope(f, c, r.get("scope"))
        s.emit(TargetField=s.target(f, sanitize(c), "Diversity"), SourceField=s.column(f),
               MatchValue=c, QuantificationType="diversity", ScopeField=sf, ScopeMatchValue=sv)

def a_activation(s, r):
    f, act = r["field"], s.cfg["activation"]
    for c in s.select(f, r["select"]):
        sf, sv = s.scope(f, c, r.get("scope"))
        s.emit(TargetField=s.target(f, sanitize(c), "Activated"),
               SourceField=act["source_field"], MatchValue=act["match_value"],
               QuantificationType="sum", ScopeField=sf, ScopeMatchValue=sv)

ARCHETYPES = {"sum": a_sum, "score": a_score, "eds": a_eds,
              "pct_positive": a_pct_positive, "diversity": a_diversity, "activation": a_activation}


def build(cfg: dict) -> str:
    apply_myeloid_coarse_derivation(cfg)
    s = Spec(cfg)
    for r in cfg["rules"]:
        ARCHETYPES[r["archetype"]](s, r)
    return "\n".join(["\t".join(COLUMNS)] + ["\t".join(row) for row in s.rows]) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument(
        "--models-dir",
        type=Path,
        default=None,
        help="Load class lists from pkl in this dir (in-memory only; does not rewrite YAML)",
    )
    args = ap.parse_args()
    cfg = load_quantify_config(args.config)
    if args.models_dir is not None:
        refresh_classes_from_models(cfg, args.models_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tsv = build(cfg)
    args.out.write_text(tsv, encoding="utf-8")
    print(f"wrote {args.out} ({tsv.count(chr(10)) - 1} rows)")


if __name__ == "__main__":
    main()
