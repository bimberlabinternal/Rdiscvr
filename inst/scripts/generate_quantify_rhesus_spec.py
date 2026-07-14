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
DELIM = "__"
IMMUNE_FIELD = "RIRA_Immune_v2.cellclass"
TNK_FIELD = "RIRA_TNK_v2.cellclass"
MYE_FINE = "RIRA_Myeloid_v3.cellclass"
MYE_COARSE = "RIRA_Myeloid_v3.coarseclass"
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
TNK_DERIVED = ["Other"]
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
HARDCODED_TNK = [
    "CD4+ T Cells",
    "CD8+ T Cells",
    "Gamma/Delta Cells",
    "NK Cells",
    "Other",
    "Unassigned",
    "Ambiguous",
    "Unknown",
]
HARDCODED_MYE_FINE = [
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
HARDCODED_MYE_COARSE = [
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

TNK_SCORE_MODULES = [
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
    return DELIM.join(parts)


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
    tnk_classes: list[str],
    myeloid_fine: list[str],
    myeloid_coarse: list[str],
) -> list[str]:
    rows = [SPEC_HEADER]

    # Block 1 — Immune (no scope): all sum rows for RIRA_Immune_v2.cellclass
    for cls in immune_classes:
        rows.append(
            spec_row(join_target("Immune", sanitize(cls)), IMMUNE_FIELD, cls, "sum")
        )

    # Block 2 — TNK (scoped T_NK): sums
    for cls in tnk_classes:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cls)),
                TNK_FIELD,
                cls,
                "sum",
                scope_field=SCOPE_IMMUNE,
                scope_match="T_NK",
            )
        )

    # Block 2a — TNK score modules (scoped T_NK)
    for score_suffix, score_field in TNK_SCORE_MODULES:
        for cls in tnk_classes:
            rows.append(
                spec_row(
                    join_target("TNK", sanitize(cls), score_suffix),
                    TNK_FIELD,
                    cls,
                    "score",
                    score_field,
                    "quantiles",
                    SCOPE_IMMUNE,
                    "T_NK",
                )
            )

    # Block 2b — CD4/CD8 EDS phenotypes (scoped T_NK)
    for cls in EDS_T_CELL_CLASSES:
        for phenotype in EDS_PHENOTYPES:
            rows.append(
                spec_row(
                    join_target("TNK", sanitize(cls), phenotype),
                    TNK_FIELD,
                    cls,
                    "sum",
                    scope_field=SCOPE_IMMUNE,
                    scope_match="T_NK",
                    effector_differentiation_score_field=EDS_METADATA_FIELD,
                    effector_differentiation_cutpoint_low=EDS_CUTPOINT_LOW,
                    effector_differentiation_cutpoint_high=EDS_CUTPOINT_HIGH,
                    subset_phenotype_output_field_name=phenotype,
                )
            )

    # Block 2c — CD4/CD8 EDS score quantiles (scoped T_NK)
    for cls in EDS_T_CELL_CLASSES:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cls), EDS_SCORE_SUFFIX),
                TNK_FIELD,
                cls,
                "score",
                EDS_METADATA_FIELD,
                "quantiles",
                SCOPE_IMMUNE,
                "T_NK",
            )
        )

    # Block 2d — CD4/CD8 activation sums (scoped by TNK class)
    for cls in ACTIVATION_T_CELL_CLASSES:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cls), "Activated"),
                ACTIVATION_SOURCE_FIELD,
                ACTIVATION_MATCH_VALUE,
                "sum",
                scope_field=TNK_FIELD,
                scope_match=cls,
            )
        )

    # Block 2e — pct_positive genes by TNK class (scoped T_NK)
    for gene in PCT_POSITIVE_GENES:
        for cls in pct_positive_classes_for_gene(gene):
            rows.append(
                spec_row(
                    join_target("TNK", sanitize(cls), "PctPositive", gene),
                    TNK_FIELD,
                    cls,
                    "pct_positive",
                    gene,
                    scope_field=SCOPE_IMMUNE,
                    scope_match="T_NK",
                )
            )

    # Block 2f — TCR diversity for CD4+ and CD8+ (scoped T_NK)
    for cls in DIVERSITY_T_CELL_CLASSES:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cls), "Diversity"),
                TNK_FIELD,
                cls,
                "diversity",
                scope_field=SCOPE_IMMUNE,
                scope_match="T_NK",
            )
        )

    # Block 3 — Myeloid (scoped Myeloid): fine sums, coarse sums, then IFN on coarse
    for cls in myeloid_fine:
        rows.append(
            spec_row(
                join_target("Myeloid", "fine", sanitize(cls)),
                MYE_FINE,
                cls,
                "sum",
                scope_field=SCOPE_IMMUNE,
                scope_match="Myeloid",
            )
        )
    for cls in myeloid_coarse:
        rows.append(
            spec_row(
                join_target("Myeloid", "coarse", sanitize(cls)),
                MYE_COARSE,
                cls,
                "sum",
                scope_field=SCOPE_IMMUNE,
                scope_match="Myeloid",
            )
        )
    for cls in myeloid_coarse:
        rows.append(
            spec_row(
                join_target(
                    "Myeloid", "coarse", sanitize(cls), INTERFERON_RESPONSE_SCORE_SUFFIX
                ),
                MYE_COARSE,
                cls,
                "score",
                "Interferon_Response_UCell",
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
            list(HARDCODED_TNK),
            list(HARDCODED_MYE_FINE),
            list(HARDCODED_MYE_COARSE),
        )

    immune_pkl = load_pkl_classes(models_dir / "RIRA_Immune_v2.pkl")
    tnk_pkl = load_pkl_classes(models_dir / "RIRA_TNK_v2.pkl")
    myeloid_pkl = load_pkl_classes(models_dir / "RIRA_FineScope_Myeloid_v3.pkl")

    immune_classes = immune_pkl + [c for c in IMMUNE_DERIVED if c not in immune_pkl]
    tnk_classes = tnk_pkl + TNK_DERIVED + [c for c in PIPELINE_META if c not in tnk_pkl]
    myeloid_fine = myeloid_pkl + [c for c in PIPELINE_META if c not in myeloid_pkl]
    myeloid_coarse = derive_coarse(myeloid_pkl)
    return immune_classes, tnk_classes, myeloid_fine, myeloid_coarse


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

    immune_classes, tnk_classes, myeloid_fine, myeloid_coarse = resolve_class_lists(
        args.models_dir, args.classes_only
    )
    rows = build_spec_rows(immune_classes, tnk_classes, myeloid_fine, myeloid_coarse)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"Immune classes ({len(immune_classes)}): {immune_classes}")
    print(f"TNK classes ({len(tnk_classes)}): {tnk_classes}")
    print(f"Myeloid fine ({len(myeloid_fine)}): {myeloid_fine}")
    print(f"Myeloid coarse ({len(myeloid_coarse)}): {myeloid_coarse}")
    print(f"Total spec rows (excl. header): {len(rows) - 1}")


if __name__ == "__main__":
    main()
