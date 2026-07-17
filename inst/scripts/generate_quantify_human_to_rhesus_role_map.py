#!/usr/bin/env python3
"""Generate inst/extdata/quantify_human_to_rhesus_role_map.tsv.

Maps CellTypist Immune_All High/Low labels to RIRA metadata fields for
coerceToRIRA quantification. Class lists and Low→High parent come from
quantify_human_immune.yml; coerce rules from quantify_human_to_rhesus_role_map.yml.

Column semantics:
  - rhesusLabel on RIRA_Immune_v2.cellclass: coarse RIRA compartment.
    T_NK is a mixed lymphoid bin (T + NK + helper ILC), not "T cells only".
  - rhesusLabel on RIRA_TNK_v2.cellclass / Myeloid fields: fine class used
    by the rhesus quantify spec after coerce.
  - lineageRole: documentation-only semantic tag (T, NK_ILC, CD4_like, …).
    .ApplyHumanToRhesusRoleMap does not read this column at runtime.

Usage (from package root):
  python3 inst/scripts/generate_quantify_human_to_rhesus_role_map.py
"""
from __future__ import annotations

from pathlib import Path

import yaml

from build_quantify_spec import load_quantify_config

SCRIPT_DIR = Path(__file__).resolve().parent
HUMAN_YAML = SCRIPT_DIR / "quantify_human_immune.yml"
ROLE_MAP_YAML = SCRIPT_DIR / "quantify_human_to_rhesus_role_map.yml"
DEFAULT_OUT = Path("inst/extdata/quantify_human_to_rhesus_role_map.tsv")


def load_role_map_config(path: Path) -> dict:
    return yaml.safe_load(path.read_text())


def _expand_low_tnk_map(cfg: dict) -> dict[str, tuple[str, str]]:
    """Role-set labels → (RIRA TNK class, lineageRole)."""
    out: dict[str, tuple[str, str]] = {}
    for role_name, (tnk_label, lineage) in cfg["tnk_role_map"].items():
        for label in cfg[role_name]:
            out[label] = (tnk_label, lineage)
    return out


def _parse_cycling_entry(entry: list) -> tuple[str, str, str | None, str | None]:
    immune, role, tnk, coarse = entry
    return immune, role, None if tnk == "null" else tnk, None if coarse == "null" else coarse


def immune_from_low_parent(
    low_label: str, parent_map: dict, high_immune_map: dict
) -> tuple[str, str] | None:
    parent = parent_map.get(low_label)
    if not parent:
        return None
    mapped = high_immune_map.get(parent)
    if not mapped:
        return None
    return mapped[0], mapped[1]


def build_rows(human_cfg: dict, role_cfg: dict) -> list[str]:
    high_classes = human_cfg["fields"]["high"]["classes"]
    low_classes = human_cfg["fields"]["low"]["classes"]
    parent_map = human_cfg["fields"]["low"]["parent"]
    high_immune_map = {
        k: (v[0], v[1]) for k, v in role_cfg["high_immune_map"].items()
    }
    high_source = role_cfg["sources"]["high"]
    low_source = role_cfg["sources"]["low"]
    targets = role_cfg["targets"]
    low_tnk_map = _expand_low_tnk_map(role_cfg)
    low_cycling_map = {
        k: _parse_cycling_entry(v) for k, v in role_cfg["low_cycling_map"].items()
    }
    helper_ilc_like = set(role_cfg["helper_ilc_like"])
    ilc_parents = set(role_cfg["ilc_parents"])
    myeloid_parents = set(role_cfg["myeloid_parents"])
    b_parents = set(role_cfg["b_parents"])
    low_myeloid_coarse = role_cfg["low_myeloid_coarse"]
    low_myeloid_fine = role_cfg["low_myeloid_fine"]

    header = "\t".join(
        ["humanSourceField", "humanLabel", "rhesusTargetField", "rhesusLabel", "lineageRole"]
    )
    rows = [header]
    seen: set[tuple[str, str, str, str]] = set()

    def add(source: str, label: str, target: str, rhesus: str, role: str) -> None:
        key = (source, label, target, rhesus)
        if key in seen:
            return
        seen.add(key)
        rows.append("\t".join([source, label, target, rhesus, role]))

    for high in high_classes:
        mapped = high_immune_map.get(high)
        if mapped:
            add(high_source, high, targets["immune"], mapped[0], mapped[1])

    for low in low_classes:
        parent = parent_map.get(low)

        if low in low_cycling_map:
            immune_label, role, tnk_label, myeloid_coarse = low_cycling_map[low]
            add(low_source, low, targets["immune"], immune_label, role)
            if tnk_label:
                add(low_source, low, targets["tnk"], tnk_label, role)
            if myeloid_coarse:
                add(low_source, low, targets["myeloid_coarse"], myeloid_coarse, role)
            continue

        if (
            parent in ({"T cells"} | ilc_parents) or low in helper_ilc_like
        ) and low in low_tnk_map:
            t_and_nk_label, role = low_tnk_map[low]
            add(low_source, low, targets["immune"], "T_NK", role)
            add(low_source, low, targets["tnk"], t_and_nk_label, role)
        elif parent in b_parents:
            add(low_source, low, targets["immune"], "Bcell", "B")
        elif parent in myeloid_parents:
            immune = immune_from_low_parent(low, parent_map, high_immune_map) or (
                "Myeloid",
                "Myeloid_APC",
            )
            add(low_source, low, targets["immune"], immune[0], immune[1])
            coarse = low_myeloid_coarse.get(low)
            if coarse:
                add(low_source, low, targets["myeloid_coarse"], coarse, "Myeloid_APC")
            fine = low_myeloid_fine.get(low)
            if fine:
                add(low_source, low, targets["myeloid_fine"], fine, "Myeloid_APC")
        else:
            immune = immune_from_low_parent(low, parent_map, high_immune_map)
            if immune:
                add(low_source, low, targets["immune"], immune[0], immune[1])

    return rows


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-o",
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help="Output TSV path (default: inst/extdata/quantify_human_to_rhesus_role_map.tsv)",
    )
    args = ap.parse_args()
    human_cfg = load_quantify_config(HUMAN_YAML)
    role_cfg = load_role_map_config(ROLE_MAP_YAML)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = build_rows(human_cfg, role_cfg)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows) - 1} map rows to {args.out}")


if __name__ == "__main__":
    main()
