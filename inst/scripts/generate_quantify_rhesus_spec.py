#!/usr/bin/env python3
"""Generate inst/extdata/quantify_rhesus_spec.tsv from RIRA CellTypist pkl models.

Class lists are read from sklearn LogisticRegression.classes_ inside each pkl
(BimberLab/RIRA inst/models). Myeloid coarseclass levels follow Classify_Myeloid
in RIRA/R/CellTypist.R. Immune Non-Immune/Unknown and TNK Other are pipeline
derived labels that appear in metadata after classification.

Usage (from package root):
  python3 inst/scripts/generate_quantify_rhesus_spec.py \\
    --models-dir /path/to/RIRA/inst/models \\
    --out inst/extdata/quantify_rhesus_spec.tsv

  python3 inst/scripts/generate_quantify_rhesus_spec.py \\
    --classes-only \\
    --out inst/extdata/quantify_rhesus_spec.tsv
"""
from __future__ import annotations

import argparse
import pickle
import re
import warnings
from pathlib import Path

GROUPING = "cDNA_ID"
DELIMITER = "__"
IMMUNE_FIELD = "RIRA_Immune_v2.cellclass"
T_AND_NK_FIELD = "RIRA_TNK_v2.cellclass"
MYELOID_FINE_FIELD = "RIRA_Myeloid_v3.cellclass"
MYELOID_COARSE_FIELD = "RIRA_Myeloid_v3.coarseclass"
SCOPE_IMMUNE = IMMUNE_FIELD
EDS_METADATA_FIELD = "Tcell_EffectorDifferentiation"
EDS_CUTPOINT_LOW = "2"
EDS_CUTPOINT_HIGH = "6"
CYTOTOXICITY_SCORE_SUFFIX = "CytotoxicityScore"
INTERFERON_RESPONSE_SCORE_SUFFIX = "InterferonResponseScore"
ACTIVATION_SOURCE_FIELD = "Is_TCR_Stimulated"
ACTIVATION_MATCH_VALUE = "TRUE"

# Post-pkl labels written by RIRA Classify_* functions (see CellTypist.R).
IMMUNE_DERIVED = ["Non-Immune", "Unknown"]
T_AND_NK_DERIVED = ["Other"]
PIPELINE_META = ["Unassigned", "Ambiguous", "Unknown"]

MONOCYTE_FINE = {"CD14+ Monocytes", "CD16+ Monocytes", "Inflammatory Monocytes"}
DC_FINE = {"DC", "Mature DC"}

HARDCODED_IMMUNE = [
    "Bcell",
    "Epithelial",
    "Erythrocyte",
    "Myeloid",
    "Stromal",
    "T_NK",
    "Non-Immune",
    "Unknown",
]
HARDCODED_T_AND_NK = [
    "CD4+ T Cells",
    "CD8+ T Cells",
    "Gamma/Delta Cells",
    "NK Cells",
    "Other",
    "Unassigned",
    "Ambiguous",
    "Unknown",
]
HARDCODED_MYELOID_FINE = [
    "Alv. mac.",
    "CD14+ Monocytes",
    "CD16+ Monocytes",
    "DC",
    "Inflammatory Monocytes",
    "MacrophageM1",
    "Macrophages",
    "Mature DC",
    "Myelocytes",
    "Promyelocytes",
    "pDC",
    "Unassigned",
    "Ambiguous",
    "Unknown",
]
HARDCODED_MYELOID_COARSE = [
    "Alv. mac.",
    "Ambiguous",
    "DC",
    "MacrophageM1",
    "Macrophages",
    "Monocytes",
    "Myelocytes",
    "pDC",
    "Promyelocytes",
    "Unassigned",
    "Unknown",
]

T_AND_NK_SCORE_MODULES = [
    (CYTOTOXICITY_SCORE_SUFFIX, "Cytotoxicity_UCell"),
    (INTERFERON_RESPONSE_SCORE_SUFFIX, "Interferon_Response_UCell"),
    ("CytokineScore", "TandNK_Activation2_UCell"),
    ("MHCIIScore", "MHCII_UCell"),
    ("ProliferationScore", "Proliferation_UCell"),
    ("PerforinScore", "Perforin_UCell"),
]

EDS_T_CELL_CLASSES = ["CD4+ T Cells", "CD8+ T Cells"]
EDS_PHENOTYPES = ["Naive", "MemoryLike", "Effector"]
EDS_SCORE_SUFFIX = "EDSScore"

ACTIVATION_T_CELL_CLASSES = ["CD4+ T Cells", "CD8+ T Cells"]

# Keep in sync with .RegisterQuantifyGeneSets() in R/Quantify.R
# (RIRA gene set Quantify.PctPositive).
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
PCT_POSITIVE_EFFECTOR_CLASSES = [
    "CD4+ T Cells",
    "CD8+ T Cells",
    "NK Cells",
    "Gamma/Delta Cells",
]

DIVERSITY_T_CELL_CLASSES = ["CD4+ T Cells", "CD8+ T Cells"]

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


def load_pkl_classes(pkl_path: Path) -> list[str]:
    with pkl_path.open("rb") as handle:
        obj = pickle.load(handle)
    return list(obj["Model"].classes_)


def derive_coarse(fine_classes: list[str]) -> list[str]:
    coarse = set()
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


def sanitize(label: str) -> str:
    s = label.replace("+", "plus")
    s = re.sub(r"[^\w]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def join_target(*parts: str) -> str:
    return DELIMITER.join(parts)


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


def pct_positive_classes_for_gene(gene: str) -> list[str]:
    if gene in PCT_POSITIVE_CD4_ONLY:
        return ["CD4+ T Cells"]
    if gene in PCT_POSITIVE_NK_ONLY:
        return ["NK Cells"]
    return PCT_POSITIVE_EFFECTOR_CLASSES


def build_spec_rows(
    immune_classes: list[str],
    t_and_nk_classes: list[str],
    myeloid_fine: list[str],
    myeloid_coarse: list[str],
) -> list[str]:
    rows = [SPEC_HEADER]

    # Block 1 - Immune (no scope): all sum rows for RIRA_Immune_v2.cellclass
    for cell_class in immune_classes:
        rows.append(
            spec_row(join_target("Immune", sanitize(cell_class)), IMMUNE_FIELD, cell_class, "sum")
        )

    # Block 2 - TNK (scoped T_NK): sums
    for cell_class in t_and_nk_classes:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cell_class)),
                T_AND_NK_FIELD,
                cell_class,
                "sum",
                scope_field=SCOPE_IMMUNE,
                scope_match="T_NK",
            )
        )

    # Block 2a - TNK score modules (scoped T_NK)
    for score_suffix, score_field in T_AND_NK_SCORE_MODULES:
        for cell_class in t_and_nk_classes:
            rows.append(
                spec_row(
                    join_target("TNK", sanitize(cell_class), score_suffix),
                    T_AND_NK_FIELD,
                    cell_class,
                    "score",
                    score_field,
                    "quantiles",
                    SCOPE_IMMUNE,
                    "T_NK",
                )
            )

    # Block 2b - CD4/CD8 EDS phenotypes (scoped T_NK)
    for cell_class in EDS_T_CELL_CLASSES:
        for phenotype in EDS_PHENOTYPES:
            rows.append(
                spec_row(
                    join_target("TNK", sanitize(cell_class), phenotype),
                    T_AND_NK_FIELD,
                    cell_class,
                    "sum",
                    scope_field=SCOPE_IMMUNE,
                    scope_match="T_NK",
                    effector_differentiation_score_field=EDS_METADATA_FIELD,
                    effector_differentiation_cutpoint_low=EDS_CUTPOINT_LOW,
                    effector_differentiation_cutpoint_high=EDS_CUTPOINT_HIGH,
                    subset_phenotype_output_field_name=phenotype,
                )
            )

    # Block 2c - CD4/CD8 EDS score quantiles (scoped T_NK)
    for cell_class in EDS_T_CELL_CLASSES:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cell_class), EDS_SCORE_SUFFIX),
                T_AND_NK_FIELD,
                cell_class,
                "score",
                EDS_METADATA_FIELD,
                "quantiles",
                SCOPE_IMMUNE,
                "T_NK",
            )
        )

    # Block 2d - CD4/CD8 activation sums (scoped by TNK class)
    for cell_class in ACTIVATION_T_CELL_CLASSES:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cell_class), "Activated"),
                ACTIVATION_SOURCE_FIELD,
                ACTIVATION_MATCH_VALUE,
                "sum",
                scope_field=T_AND_NK_FIELD,
                scope_match=cell_class,
            )
        )

    # Block 2e - pct_positive genes by TNK class (scoped T_NK)
    for gene in PCT_POSITIVE_GENES:
        for cell_class in pct_positive_classes_for_gene(gene):
            rows.append(
                spec_row(
                    join_target("TNK", sanitize(cell_class), "PctPositive", gene),
                    T_AND_NK_FIELD,
                    cell_class,
                    "pct_positive",
                    gene,
                    scope_field=SCOPE_IMMUNE,
                    scope_match="T_NK",
                )
            )

    # Block 2f - TCR diversity for CD4+ and CD8+ (scoped T_NK)
    for cell_class in DIVERSITY_T_CELL_CLASSES:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cell_class), "Diversity"),
                T_AND_NK_FIELD,
                cell_class,
                "diversity",
                scope_field=SCOPE_IMMUNE,
                scope_match="T_NK",
            )
        )

    # Block 1a - Immune parent scores (T_NK full module set; Myeloid/B IFN+MHCII)
    immune_parent_modules = {
        "T_NK": T_AND_NK_SCORE_MODULES,
        "Myeloid": [
            (INTERFERON_RESPONSE_SCORE_SUFFIX, "Interferon_Response_UCell"),
            ("MHCIIScore", "MHCII_UCell"),
        ],
        "Bcell": [
            (INTERFERON_RESPONSE_SCORE_SUFFIX, "Interferon_Response_UCell"),
            ("MHCIIScore", "MHCII_UCell"),
            ("ProliferationScore", "Proliferation_UCell"),
        ],
    }
    for immune_class, modules in immune_parent_modules.items():
        if immune_class not in immune_classes:
            continue
        for score_suffix, score_field in modules:
            rows.append(
                spec_row(
                    join_target("Immune", sanitize(immune_class), score_suffix),
                    IMMUNE_FIELD,
                    immune_class,
                    "score",
                    score_field,
                    "quantiles",
                )
            )

    # Block 3 - Myeloid (scoped Myeloid): fine sums, coarse sums, then IFN+MHCII on coarse
    for cell_class in myeloid_fine:
        rows.append(
            spec_row(
                join_target("Myeloid", "fine", sanitize(cell_class)),
                MYELOID_FINE_FIELD,
                cell_class,
                "sum",
                scope_field=SCOPE_IMMUNE,
                scope_match="Myeloid",
            )
        )
    for cell_class in myeloid_coarse:
        rows.append(
            spec_row(
                join_target("Myeloid", "coarse", sanitize(cell_class)),
                MYELOID_COARSE_FIELD,
                cell_class,
                "sum",
                scope_field=SCOPE_IMMUNE,
                scope_match="Myeloid",
            )
        )
    for cell_class in myeloid_coarse:
        for score_suffix, score_field in [
            (INTERFERON_RESPONSE_SCORE_SUFFIX, "Interferon_Response_UCell"),
            ("MHCIIScore", "MHCII_UCell"),
        ]:
            rows.append(
                spec_row(
                    join_target("Myeloid", "coarse", sanitize(cell_class), score_suffix),
                    MYELOID_COARSE_FIELD,
                    cell_class,
                    "score",
                    score_field,
                    "quantiles",
                    SCOPE_IMMUNE,
                    "Myeloid",
                )
            )

    return rows


def resolve_class_lists(models_dir: Path, classes_only: bool) -> tuple[
    list[str], list[str], list[str], list[str]
]:
    if classes_only:
        return (
            list(HARDCODED_IMMUNE),
            list(HARDCODED_T_AND_NK),
            list(HARDCODED_MYELOID_FINE),
            list(HARDCODED_MYELOID_COARSE),
        )

    immune_pkl = load_pkl_classes(models_dir / "RIRA_Immune_v2.pkl")
    t_and_nk_pkl = load_pkl_classes(models_dir / "RIRA_TNK_v2.pkl")
    myeloid_pkl = load_pkl_classes(models_dir / "RIRA_FineScope_Myeloid_v3.pkl")

    immune_classes = immune_pkl + [c for c in IMMUNE_DERIVED if c not in immune_pkl]
    t_and_nk_classes = t_and_nk_pkl + T_AND_NK_DERIVED + [c for c in PIPELINE_META if c not in t_and_nk_pkl]
    myeloid_fine = myeloid_pkl + [c for c in PIPELINE_META if c not in myeloid_pkl]
    myeloid_coarse = derive_coarse(myeloid_pkl)
    return immune_classes, t_and_nk_classes, myeloid_fine, myeloid_coarse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--models-dir",
        type=Path,
        default=Path("inst/models"),
        help="Directory containing RIRA *.pkl models",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("inst/extdata/quantify_rhesus_spec.tsv"),
    )
    parser.add_argument(
        "--classes-only",
        action="store_true",
        help="Skip pkl loading and use hardcoded class lists from the current TSV",
    )
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=UserWarning)

    immune_classes, t_and_nk_classes, myeloid_fine, myeloid_coarse = resolve_class_lists(
        args.models_dir, args.classes_only
    )
    rows = build_spec_rows(immune_classes, t_and_nk_classes, myeloid_fine, myeloid_coarse)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"Immune classes ({len(immune_classes)}): {immune_classes}")
    print(f"T/NK classes ({len(t_and_nk_classes)}): {t_and_nk_classes}")
    print(f"Myeloid fine ({len(myeloid_fine)}): {myeloid_fine}")
    print(f"Myeloid coarse ({len(myeloid_coarse)}): {myeloid_coarse}")
    print(f"Total spec rows (excl. header): {len(rows) - 1}")


if __name__ == "__main__":
    main()
