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
DELIMITER = "__"
HIGH_FIELD = "celltypist.Immune_All_High.cellclass"
LOW_FIELD = "celltypist.Immune_All_Low.cellclass"
HIGH_PKL = "Immune_All_High.pkl"
LOW_PKL = "Immune_All_Low.pkl"

PIPELINE_META = ["Unassigned", "Ambiguous", "Unknown"]

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

# Keep in sync with .RegisterQuantifyGeneSets() in R/Quantify.R
PCT_POSITIVE_GENES = [
    "PDCD1",
    "KLRK1",
    "KLRB1",
    "HAVCR2",
    "TIGIT",
    "PRF1",
    "GZMB",
    "KLRC1",
    "FCGR3",
    "FCGR3A",
    "FOXP3",
    "IL2RA",
]

PCT_POSITIVE_CD4_ONLY = {"FOXP3", "IL2RA"}
PCT_POSITIVE_NK_ONLY = {"FCGR3", "FCGR3A"}

T_SCORE_MODULES = [
    ("CytotoxicityScore", "Cytotoxicity_UCell"),
    ("InterferonResponseScore", "Interferon_Response_UCell"),
    ("CytokineScore", "TandNK_Activation2_UCell"),
    ("MHCIIScore", "MHCII_UCell"),
    ("ProliferationScore", "Proliferation_UCell"),
    ("PerforinScore", "Perforin_UCell"),
]

NK_ILC_SCORE_MODULES = [
    ("CytotoxicityScore", "Cytotoxicity_UCell"),
    ("InterferonResponseScore", "Interferon_Response_UCell"),
    ("MHCIIScore", "MHCII_UCell"),
    ("ProliferationScore", "Proliferation_UCell"),
    ("PerforinScore", "Perforin_UCell"),
]

ANTIGEN_PRESENTING_SCORE_MODULES = [
    ("InterferonResponseScore", "Interferon_Response_UCell"),
    ("MHCIIScore", "MHCII_UCell"),
    ("ProliferationScore", "Proliferation_UCell"),
]

HIGH_ANTIGEN_PRESENTING_PARENTS = {
    "B cells",
    "B-cell lineage",
    "Plasma cells",
    "Macrophages",
    "Monocytes",
    "DC",
    "MNP",
    "Mono-mac",
    "Myelocytes",
    "Promyelocytes",
    "pDC",
    "Mast cells",
    "Granulocytes",
}

LOW_ANTIGEN_PRESENTING_PARENTS = {
    "B cells",
    "B-cell lineage",
    "Plasma cells",
    "Macrophages",
    "Monocytes",
    "DC",
    "MNP",
    "Mono-mac",
    "Myelocytes",
    "Promyelocytes",
    "pDC",
    "Mast cells",
    "Granulocytes",
}

CD4_LIKE = {
    "Tcm/Naive helper T cells",
    "Tem/Effector helper T cells",
    "Tem/Effector helper T cells PD1+",
    "Memory CD4+ cytotoxic T cells",
    "Regulatory T cells",
    "Follicular helper T cells",
    "Type 1 helper T cells",
    "Type 17 helper T cells",
    "Treg(diff)",
    "T(agonist)",
}

CD8_LIKE = {
    "Tcm/Naive cytotoxic T cells",
    "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells",
    "Trm cytotoxic T cells",
    "CD8a/a",
    "CD8a/b(entry)",
    "MAIT cells",
}

EDS_PHENOTYPES = ["Naive", "MemoryLike", "Effector"]
EDS_METADATA_FIELD = "Tcell_EffectorDifferentiation"
EDS_CUTPOINT_LOW = "2"
EDS_CUTPOINT_HIGH = "6"
EDS_SCORE_SUFFIX = "EDSScore"

ACTIVATION_SOURCE_FIELD = "Is_TCR_Stimulated"
ACTIVATION_MATCH_VALUE = "TRUE"

PCT_POSITIVE_EFFECTOR_CLASSES = [
    "Tcm/Naive cytotoxic T cells",
    "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells",
    "Trm cytotoxic T cells",
    "CD8a/a",
    "CD8a/b(entry)",
    "MAIT cells",
    "Tcm/Naive helper T cells",
    "Tem/Effector helper T cells",
    "Tem/Effector helper T cells PD1+",
    "Memory CD4+ cytotoxic T cells",
    "Regulatory T cells",
    "Follicular helper T cells",
    "Type 1 helper T cells",
    "Type 17 helper T cells",
    "Treg(diff)",
    "T(agonist)",
    "gamma-delta T cells",
    "CRTAM+ gamma-delta T cells",
    "Cycling gamma-delta T cells",
    "NK cells",
    "CD16+ NK cells",
    "CD16- NK cells",
    "Cycling NK cells",
    "Transitional NK",
    "NKT cells",
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
    return DELIMITER.join(parts)


def load_pkl_classes(pkl_path: Path) -> list[str]:
    with pkl_path.open("rb") as handle:
        obj = pickle.load(handle)
    return list(obj["Model"].classes_)


def with_pipeline_meta(classes: list[str]) -> list[str]:
    return classes + [c for c in PIPELINE_META if c not in classes]


def lows_with_parent(low_classes: list[str], parent: str) -> list[str]:
    return [c for c in low_classes if LOW_HIGH_SCOPE.get(c) == parent]


def pct_positive_classes_for_gene(gene: str) -> list[str]:
    if gene in PCT_POSITIVE_CD4_ONLY:
        return [c for c in CD4_LIKE if c in PCT_POSITIVE_EFFECTOR_CLASSES]
    if gene in PCT_POSITIVE_NK_ONLY:
        return [
            c
            for c in PCT_POSITIVE_EFFECTOR_CLASSES
            if c
            in {
                "NK cells",
                "CD16+ NK cells",
                "CD16- NK cells",
                "Cycling NK cells",
                "Transitional NK",
            }
        ]
    return [c for c in PCT_POSITIVE_EFFECTOR_CLASSES if c in CD4_LIKE | CD8_LIKE | {
        "gamma-delta T cells",
        "CRTAM+ gamma-delta T cells",
        "Cycling gamma-delta T cells",
        "NK cells",
        "CD16+ NK cells",
        "CD16- NK cells",
        "Cycling NK cells",
        "Transitional NK",
        "NKT cells",
    }]


def spec_row(
    target: str,
    source: str,
    match: str,
    qtype: str,
    qsource: str = "",
    qscore: str = "",
    scope_field: str = "",
    scope_match: str = "",
    effector_differentiation_score_field: str = "",
    effector_differentiation_cutpoint_low: str = "",
    effector_differentiation_cutpoint_high: str = "",
    subset_phenotype_output_field_name: str = "",
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
            effector_differentiation_score_field,
            effector_differentiation_cutpoint_low,
            effector_differentiation_cutpoint_high,
            subset_phenotype_output_field_name,
        ]
    )


def append_score_rows(
    rows: list[str],
    prefix: str,
    source: str,
    match: str,
    modules: list[tuple[str, str]],
    scope_field: str = "",
    scope_match: str = "",
) -> None:
    for score_suffix, score_field in modules:
        rows.append(
            spec_row(
                join_target(prefix, sanitize(match), score_suffix),
                source,
                match,
                "score",
                score_field,
                "quantiles",
                scope_field=scope_field,
                scope_match=scope_match,
            )
        )


def build_spec_rows(high_classes: list[str], low_classes: list[str]) -> list[str]:
    rows = [SPEC_HEADER]
    t_cell_labels = lows_with_parent(low_classes, "T cells")
    ilc_labels = lows_with_parent(low_classes, "ILC")
    cd4_labels = [c for c in t_cell_labels if c in CD4_LIKE]
    cd8_labels = [c for c in t_cell_labels if c in CD8_LIKE]
    diversity_labels = cd4_labels + cd8_labels

    # Block 1 - Immune High sums
    for cell_class in high_classes:
        rows.append(
            spec_row(join_target("ImmuneHigh", sanitize(cell_class)), HIGH_FIELD, cell_class, "sum")
        )

    # Block 2 - Immune Low sums
    for cell_class in low_classes:
        scope_match = LOW_HIGH_SCOPE.get(cell_class, "")
        scope_field = HIGH_FIELD if scope_match else ""
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cell_class)),
                LOW_FIELD,
                cell_class,
                "sum",
                scope_field=scope_field,
                scope_match=scope_match,
            )
        )

    # Block 3 - High parent scores (T, ILC, APC parents)
    for cell_class in high_classes:
        if cell_class == "T cells":
            append_score_rows(rows, "ImmuneHigh", HIGH_FIELD, cell_class, T_SCORE_MODULES)
        elif cell_class == "ILC":
            append_score_rows(rows, "ImmuneHigh", HIGH_FIELD, cell_class, NK_ILC_SCORE_MODULES)
        elif cell_class in HIGH_ANTIGEN_PRESENTING_PARENTS:
            append_score_rows(rows, "ImmuneHigh", HIGH_FIELD, cell_class, ANTIGEN_PRESENTING_SCORE_MODULES)

    # Block 4 - Low T-cell score modules
    for cell_class in t_cell_labels:
        append_score_rows(
            rows,
            "ImmuneLow",
            LOW_FIELD,
            cell_class,
            T_SCORE_MODULES,
            scope_field=HIGH_FIELD,
            scope_match="T cells",
        )

    # Block 5 - Low ILC/NK score modules
    for cell_class in ilc_labels:
        append_score_rows(
            rows,
            "ImmuneLow",
            LOW_FIELD,
            cell_class,
            NK_ILC_SCORE_MODULES,
            scope_field=HIGH_FIELD,
            scope_match="ILC",
        )

    # Block 6 - Low APC score modules
    for cell_class in low_classes:
        parent = LOW_HIGH_SCOPE.get(cell_class, "")
        if parent in LOW_ANTIGEN_PRESENTING_PARENTS:
            append_score_rows(
                rows,
                "ImmuneLow",
                LOW_FIELD,
                cell_class,
                ANTIGEN_PRESENTING_SCORE_MODULES,
                scope_field=HIGH_FIELD,
                scope_match=parent,
            )

    # Block 7 - EDS phenotypes + score on CD4/CD8-like Low labels
    for cell_class in cd4_labels + cd8_labels:
        for phenotype in EDS_PHENOTYPES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cell_class), phenotype),
                    LOW_FIELD,
                    cell_class,
                    "sum",
                    scope_field=HIGH_FIELD,
                    scope_match="T cells",
                    effector_differentiation_score_field=EDS_METADATA_FIELD,
                    effector_differentiation_cutpoint_low=EDS_CUTPOINT_LOW,
                    effector_differentiation_cutpoint_high=EDS_CUTPOINT_HIGH,
                    subset_phenotype_output_field_name=phenotype,
                )
            )
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cell_class), EDS_SCORE_SUFFIX),
                LOW_FIELD,
                cell_class,
                "score",
                EDS_METADATA_FIELD,
                "quantiles",
                scope_field=HIGH_FIELD,
                scope_match="T cells",
            )
        )

    # Block 8 - pct_positive genes with gene×role filters
    for gene in PCT_POSITIVE_GENES:
        for cell_class in pct_positive_classes_for_gene(gene):
            if cell_class not in low_classes:
                continue
            parent = LOW_HIGH_SCOPE.get(cell_class, "")
            scope_match = "T cells" if parent == "T cells" else "ILC" if parent == "ILC" else ""
            scope_field = HIGH_FIELD if scope_match else ""
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cell_class), "PctPositive", gene),
                    LOW_FIELD,
                    cell_class,
                    "pct_positive",
                    gene,
                    scope_field=scope_field,
                    scope_match=scope_match,
                )
            )

    # Block 9 - TCR diversity on all CD4/CD8-like Low labels
    for cell_class in diversity_labels:
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cell_class), "Diversity"),
                LOW_FIELD,
                cell_class,
                "diversity",
                scope_field=HIGH_FIELD,
                scope_match="T cells",
            )
        )

    # Block 10 - activation sums for Low T labels + High T parent
    rows.append(
        spec_row(
            join_target("ImmuneHigh", sanitize("T cells"), "Activated"),
            ACTIVATION_SOURCE_FIELD,
            ACTIVATION_MATCH_VALUE,
            "sum",
            scope_field=HIGH_FIELD,
            scope_match="T cells",
        )
    )
    for cell_class in t_cell_labels:
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cell_class), "Activated"),
                ACTIVATION_SOURCE_FIELD,
                ACTIVATION_MATCH_VALUE,
                "sum",
                scope_field=LOW_FIELD,
                scope_match=cell_class,
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
