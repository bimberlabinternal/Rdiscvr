#!/usr/bin/env python3
"""Generate inst/extdata/quantify_human_immune_spec.tsv from CellTypist pkl models.

Class lists are read from sklearn LogisticRegression.classes_ inside
Immune_All_High.pkl and Immune_All_Low.pkl (CellTypist Pan Immune Atlas v2).
Low→High ScopeField parents follow Basic_celltype_information.xlsx from
Teichlab/celltypist_wiki. Pipeline Unassigned/Ambiguous/Unknown labels are
appended after pkl classes (same pattern as the rhesus generator).

Usage (from package root):
  python3 inst/scripts/generate_quantify_human_spec.py \\
    --models-dir /path/to/celltypist/models \\
    --out inst/extdata/quantify_human_immune_spec.tsv

  python3 inst/scripts/generate_quantify_human_spec.py \\
    --classes-only \\
    --out inst/extdata/quantify_human_immune_spec.tsv
"""
from __future__ import annotations

import argparse
import pickle
import re
import warnings
from pathlib import Path

GROUPING = "cDNA_ID"
DELIM = "__"
HIGH_FIELD = "celltypist.Immune_All_High.cellclass"
LOW_FIELD = "celltypist.Immune_All_Low.cellclass"
HIGH_PKL = "Immune_All_High.pkl"
LOW_PKL = "Immune_All_Low.pkl"

PIPELINE_META = ["Unassigned", "Ambiguous", "Unknown"]

# Frozen snapshot of Model.classes_ from CellTypist Pan Immune Atlas v2
# (https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/).
HARDCODED_HIGH = [
    "B cells",
    "B-cell lineage",
    "Cycling cells",
    "DC",
    "DC precursor",
    "Double-negative thymocytes",
    "Double-positive thymocytes",
    "ETP",
    "Early MK",
    "Endothelial cells",
    "Epithelial cells",
    "Erythrocytes",
    "Erythroid",
    "Fibroblasts",
    "Granulocytes",
    "HSC/MPP",
    "ILC",
    "ILC precursor",
    "MNP",
    "Macrophages",
    "Mast cells",
    "Megakaryocyte precursor",
    "Megakaryocytes/platelets",
    "Mono-mac",
    "Monocyte precursor",
    "Monocytes",
    "Myelocytes",
    "Plasma cells",
    "Promyelocytes",
    "T cells",
    "pDC",
    "pDC precursor",
]

HARDCODED_LOW = [
    "Age-associated B cells",
    "Alveolar macrophages",
    "B cells",
    "CD16+ NK cells",
    "CD16- NK cells",
    "CD8a/a",
    "CD8a/b(entry)",
    "CMP",
    "CRTAM+ gamma-delta T cells",
    "Classical monocytes",
    "Cycling B cells",
    "Cycling DCs",
    "Cycling NK cells",
    "Cycling T cells",
    "Cycling gamma-delta T cells",
    "Cycling monocytes",
    "DC",
    "DC precursor",
    "DC1",
    "DC2",
    "DC3",
    "Double-negative thymocytes",
    "Double-positive thymocytes",
    "ELP",
    "ETP",
    "Early MK",
    "Early erythroid",
    "Early lymphoid/T lymphoid",
    "Endothelial cells",
    "Epithelial cells",
    "Erythrocytes",
    "Erythrophagocytic macrophages",
    "Fibroblasts",
    "Follicular B cells",
    "Follicular helper T cells",
    "GMP",
    "Germinal center B cells",
    "Granulocytes",
    "HSC/MPP",
    "Hofbauer cells",
    "ILC",
    "ILC precursor",
    "ILC1",
    "ILC2",
    "ILC3",
    "Intermediate macrophages",
    "Intestinal macrophages",
    "Kidney-resident macrophages",
    "Kupffer cells",
    "Large pre-B cells",
    "Late erythroid",
    "MAIT cells",
    "MEMP",
    "MNP",
    "Macrophages",
    "Mast cells",
    "Megakaryocyte precursor",
    "Megakaryocyte-erythroid-mast cell progenitor",
    "Megakaryocytes/platelets",
    "Memory B cells",
    "Memory CD4+ cytotoxic T cells",
    "Mid erythroid",
    "Migratory DCs",
    "Mono-mac",
    "Monocyte precursor",
    "Monocytes",
    "Myelocytes",
    "NK cells",
    "NKT cells",
    "Naive B cells",
    "Neutrophil-myeloid progenitor",
    "Neutrophils",
    "Non-classical monocytes",
    "Plasma cells",
    "Plasmablasts",
    "Pre-pro-B cells",
    "Pro-B cells",
    "Proliferative germinal center B cells",
    "Promyelocytes",
    "Regulatory T cells",
    "Small pre-B cells",
    "T(agonist)",
    "Tcm/Naive cytotoxic T cells",
    "Tcm/Naive helper T cells",
    "Tem/Effector helper T cells",
    "Tem/Effector helper T cells PD1+",
    "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells",
    "Transitional B cells",
    "Transitional DC",
    "Transitional NK",
    "Treg(diff)",
    "Trm cytotoxic T cells",
    "Type 1 helper T cells",
    "Type 17 helper T cells",
    "gamma-delta T cells",
    "pDC",
    "pDC precursor",
]

# Low→High parents from celltypist_wiki Basic_celltype_information.xlsx (v2).
LOW_HIGH_SCOPE = {
    "Age-associated B cells": "B cells",
    "Alveolar macrophages": "Macrophages",
    "B cells": "B cells",
    "CD16+ NK cells": "ILC",
    "CD16- NK cells": "ILC",
    "CD8a/a": "T cells",
    "CD8a/b(entry)": "T cells",
    "CMP": "HSC/MPP",
    "CRTAM+ gamma-delta T cells": "T cells",
    "Classical monocytes": "Monocytes",
    "Cycling B cells": "Cycling cells",
    "Cycling DCs": "Cycling cells",
    "Cycling NK cells": "Cycling cells",
    "Cycling T cells": "Cycling cells",
    "Cycling gamma-delta T cells": "Cycling cells",
    "Cycling monocytes": "Cycling cells",
    "DC": "DC",
    "DC precursor": "DC precursor",
    "DC1": "DC",
    "DC2": "DC",
    "DC3": "DC",
    "Double-negative thymocytes": "Double-negative thymocytes",
    "Double-positive thymocytes": "Double-positive thymocytes",
    "ELP": "HSC/MPP",
    "ETP": "ETP",
    "Early MK": "Early MK",
    "Early erythroid": "Erythroid",
    "Early lymphoid/T lymphoid": "HSC/MPP",
    "Endothelial cells": "Endothelial cells",
    "Epithelial cells": "Epithelial cells",
    "Erythrocytes": "Erythrocytes",
    "Erythrophagocytic macrophages": "Macrophages",
    "Fibroblasts": "Fibroblasts",
    "Follicular B cells": "B cells",
    "Follicular helper T cells": "T cells",
    "GMP": "HSC/MPP",
    "Germinal center B cells": "B cells",
    "Granulocytes": "Granulocytes",
    "HSC/MPP": "HSC/MPP",
    "Hofbauer cells": "Macrophages",
    "ILC": "ILC",
    "ILC precursor": "ILC precursor",
    "ILC1": "ILC",
    "ILC2": "ILC",
    "ILC3": "ILC",
    "Intermediate macrophages": "Macrophages",
    "Intestinal macrophages": "Macrophages",
    "Kidney-resident macrophages": "Macrophages",
    "Kupffer cells": "Macrophages",
    "Large pre-B cells": "B-cell lineage",
    "Late erythroid": "Erythroid",
    "MAIT cells": "T cells",
    "MEMP": "HSC/MPP",
    "MNP": "MNP",
    "Macrophages": "Macrophages",
    "Mast cells": "Mast cells",
    "Megakaryocyte precursor": "Megakaryocyte precursor",
    "Megakaryocyte-erythroid-mast cell progenitor": "HSC/MPP",
    "Megakaryocytes/platelets": "Megakaryocytes/platelets",
    "Memory B cells": "B cells",
    "Memory CD4+ cytotoxic T cells": "T cells",
    "Mid erythroid": "Erythroid",
    "Migratory DCs": "DC",
    "Mono-mac": "Mono-mac",
    "Monocyte precursor": "Monocyte precursor",
    "Monocytes": "Monocytes",
    "Myelocytes": "Myelocytes",
    "NK cells": "ILC",
    "NKT cells": "T cells",
    "Naive B cells": "B cells",
    "Neutrophil-myeloid progenitor": "HSC/MPP",
    "Neutrophils": "Granulocytes",
    "Non-classical monocytes": "Monocytes",
    "Plasma cells": "Plasma cells",
    "Plasmablasts": "Plasma cells",
    "Pre-pro-B cells": "B-cell lineage",
    "Pro-B cells": "B-cell lineage",
    "Proliferative germinal center B cells": "B cells",
    "Promyelocytes": "Promyelocytes",
    "Regulatory T cells": "T cells",
    "Small pre-B cells": "B-cell lineage",
    "T(agonist)": "T cells",
    "Tcm/Naive cytotoxic T cells": "T cells",
    "Tcm/Naive helper T cells": "T cells",
    "Tem/Effector helper T cells": "T cells",
    "Tem/Effector helper T cells PD1+": "T cells",
    "Tem/Temra cytotoxic T cells": "T cells",
    "Tem/Trm cytotoxic T cells": "T cells",
    "Transitional B cells": "B cells",
    "Transitional DC": "DC",
    "Transitional NK": "ILC",
    "Treg(diff)": "T cells",
    "Trm cytotoxic T cells": "T cells",
    "Type 1 helper T cells": "T cells",
    "Type 17 helper T cells": "T cells",
    "gamma-delta T cells": "T cells",
    "pDC": "pDC",
    "pDC precursor": "pDC precursor",
}

LOW_SCORE_MODULES = [
    ("CytotoxicityScore", "Cytotoxicity_UCell"),
    ("InterferonResponseScore", "Interferon_Response_UCell"),
]

LOW_PCT_POSITIVE_GENES = ["PDCD1", "KLRK1", "GZMB"]  # subset of Quantify.PctPositive in R/Quantify.R

# Curated diversity targets (must remain in the Low class list).
DIVERSITY_LOW_CLASSES = [
    "Tcm/Naive helper T cells",
    "Tem/Temra cytotoxic T cells",
]

SPEC_HEADER = "\t".join(
    [
        "GroupingVariable",
        "TargetField",
        "SourceField",
        "MatchValue",
        "QuantificationType",
        "QuantificationSourceField",
        "QuantificationScoreType",
        "ScopeField",
        "ScopeMatchValue",
        "EffectorDifferentiationScoreField",
        "EffectorDifferentiationCutpointLow",
        "EffectorDifferentiationCutpointHigh",
        "SubsetPhenotypeOutputFieldName",
    ]
)


def sanitize(label: str) -> str:
    s = label.replace("+", "plus")
    s = re.sub(r"[^\w]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def join_target(*parts: str) -> str:
    return DELIM.join(parts)


def load_pkl_classes(pkl_path: Path) -> list[str]:
    with pkl_path.open("rb") as handle:
        obj = pickle.load(handle)
    return list(obj["Model"].classes_)


def with_pipeline_meta(classes: list[str]) -> list[str]:
    return classes + [c for c in PIPELINE_META if c not in classes]


def lows_with_parent(low_classes: list[str], parent: str) -> list[str]:
    return [c for c in low_classes if LOW_HIGH_SCOPE.get(c) == parent]


def spec_row(
    target: str,
    source: str,
    match: str,
    qtype: str,
    qsource: str = "",
    qscore: str = "",
    scope_field: str = "",
    scope_match: str = "",
) -> str:
    return "\t".join(
        [
            GROUPING,
            target,
            source,
            match,
            qtype,
            qsource,
            qscore,
            scope_field,
            scope_match,
            "",
            "",
            "",
            "",
        ]
    )


def build_spec_rows(high_classes: list[str], low_classes: list[str]) -> list[str]:
    rows = [SPEC_HEADER]
    t_cell_labels = lows_with_parent(low_classes, "T cells")
    ilc_labels = lows_with_parent(low_classes, "ILC")
    diversity_labels = [c for c in DIVERSITY_LOW_CLASSES if c in low_classes]

    # Block 1 — Immune High sums
    for cls in high_classes:
        rows.append(
            spec_row(join_target("ImmuneHigh", sanitize(cls)), HIGH_FIELD, cls, "sum")
        )

    # Block 2 — Immune Low sums (scoped to High parent when mapped)
    for cls in low_classes:
        scope_match = LOW_HIGH_SCOPE.get(cls, "")
        scope_field = HIGH_FIELD if scope_match else ""
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cls)),
                LOW_FIELD,
                cls,
                "sum",
                scope_field=scope_field,
                scope_match=scope_match,
            )
        )

    # Block 3 — T-cell score modules (all Low labels under High "T cells")
    for cls in t_cell_labels:
        for score_suffix, score_field in LOW_SCORE_MODULES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cls), score_suffix),
                    LOW_FIELD,
                    cls,
                    "score",
                    score_field,
                    "quantiles",
                    scope_field=HIGH_FIELD,
                    scope_match="T cells",
                )
            )

    # Block 4 — ILC/NK score modules (all Low labels under High "ILC")
    for cls in ilc_labels:
        for score_suffix, score_field in LOW_SCORE_MODULES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cls), score_suffix),
                    LOW_FIELD,
                    cls,
                    "score",
                    score_field,
                    "quantiles",
                    scope_field=HIGH_FIELD,
                    scope_match="ILC",
                )
            )

    # Block 5 — pct_positive on T-cell Low labels
    for cls in t_cell_labels:
        for gene in LOW_PCT_POSITIVE_GENES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cls), "PctPositive", gene),
                    LOW_FIELD,
                    cls,
                    "pct_positive",
                    gene,
                    scope_field=HIGH_FIELD,
                    scope_match="T cells",
                )
            )

    # Block 6 — TCR diversity on selected Low T-cell labels
    for cls in diversity_labels:
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cls), "Diversity"),
                LOW_FIELD,
                cls,
                "diversity",
                scope_field=HIGH_FIELD,
                scope_match="T cells",
            )
        )

    return rows


def resolve_class_lists(models_dir: Path, classes_only: bool) -> tuple[list[str], list[str]]:
    if classes_only:
        return with_pipeline_meta(list(HARDCODED_HIGH)), with_pipeline_meta(list(HARDCODED_LOW))

    high_pkl = load_pkl_classes(models_dir / HIGH_PKL)
    low_pkl = load_pkl_classes(models_dir / LOW_PKL)
    return with_pipeline_meta(high_pkl), with_pipeline_meta(low_pkl)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--models-dir",
        type=Path,
        default=Path("inst/models"),
        help="Directory containing Immune_All_High.pkl and Immune_All_Low.pkl",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("inst/extdata/quantify_human_immune_spec.tsv"),
    )
    parser.add_argument(
        "--classes-only",
        action="store_true",
        help="Skip pkl loading and use hardcoded Pan Immune Atlas v2 class lists",
    )
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=UserWarning)

    high_classes, low_classes = resolve_class_lists(args.models_dir, args.classes_only)
    unscoped = [c for c in low_classes if c not in PIPELINE_META and c not in LOW_HIGH_SCOPE]
    if unscoped:
        raise SystemExit(
            "Low labels missing from LOW_HIGH_SCOPE: " + ", ".join(unscoped)
        )

    rows = build_spec_rows(high_classes, low_classes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"High classes ({len(high_classes)}): {high_classes}")
    print(f"Low classes ({len(low_classes)}): {low_classes}")
    print(f"Total spec rows (excl. header): {len(rows) - 1}")


if __name__ == "__main__":
    main()
