#!/usr/bin/env python3
"""Generate inst/extdata/quantify_default_spec.tsv from RIRA CellTypist pkl models.

Class lists are read from sklearn LogisticRegression.classes_ inside each pkl
(BimberLab/RIRA inst/models). Myeloid coarseclass levels follow Classify_Myeloid
in RIRA/R/CellTypist.R. Immune Non-Immune/Unknown and TNK Other are pipeline
derived labels that appear in metadata after classification.

Usage (from package root):
  python3 inst/scripts/generate_quantify_default_spec.py \\
    --models-dir /path/to/RIRA/inst/models \\
    --out inst/extdata/quantify_default_spec.tsv
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
EDS_T_CELL_CLASSES = ["CD4+ T Cells", "CD8+ T Cells"]
EDS_PHENOTYPES = ["Naive", "MemoryLike", "Effector"]

# Post-pkl labels written by RIRA Classify_* functions (see CellTypist.R).
IMMUNE_DERIVED = ["Non-Immune", "Unknown"]
TNK_DERIVED = ["Other"]
PIPELINE_META = ["Unassigned", "Ambiguous", "Unknown"]

MONOCYTE_FINE = {"CD14+ Monocytes", "CD16+ Monocytes", "Inflammatory Monocytes"}
DC_FINE = {"DC", "Mature DC"}


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
        default=Path("inst/extdata/quantify_default_spec.tsv"),
    )
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=UserWarning)

    immune_pkl = load_pkl_classes(args.models_dir / "RIRA_Immune_v2.pkl")
    tnk_pkl = load_pkl_classes(args.models_dir / "RIRA_TNK_v2.pkl")
    myeloid_pkl = load_pkl_classes(args.models_dir / "RIRA_FineScope_Myeloid_v3.pkl")

    immune_classes = immune_pkl + [c for c in IMMUNE_DERIVED if c not in immune_pkl]
    tnk_classes = tnk_pkl + TNK_DERIVED + [c for c in PIPELINE_META if c not in tnk_pkl]
    myeloid_fine = myeloid_pkl + [c for c in PIPELINE_META if c not in myeloid_pkl]
    myeloid_coarse = derive_coarse(myeloid_pkl)

    rows = [
        "\t".join(
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
    ]

    # Block 1 — Immune (no scope): all sum rows for RIRA_Immune_v2.cellclass
    for cls in immune_classes:
        rows.append(
            spec_row(join_target("Immune", sanitize(cls)), IMMUNE_FIELD, cls, "sum")
        )

    # Block 2 — TNK (scoped T_NK): sums, then Cytotoxicity, then Interferon
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
    for cls in tnk_classes:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cls), CYTOTOXICITY_SCORE_SUFFIX),
                TNK_FIELD,
                cls,
                "score",
                "Cytotoxicity_UCell",
                "quantiles",
                SCOPE_IMMUNE,
                "T_NK",
            )
        )
    for cls in tnk_classes:
        rows.append(
            spec_row(
                join_target("TNK", sanitize(cls), INTERFERON_RESPONSE_SCORE_SUFFIX),
                TNK_FIELD,
                cls,
                "score",
                "Interferon_Response_UCell",
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

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"Immune classes ({len(immune_classes)}): {immune_classes}")
    print(f"TNK classes ({len(tnk_classes)}): {tnk_classes}")
    print(f"Myeloid fine ({len(myeloid_fine)}): {myeloid_fine}")
    print(f"Myeloid coarse ({len(myeloid_coarse)}): {myeloid_coarse}")
    print(f"Total spec rows (excl. header): {len(rows) - 1}")


if __name__ == "__main__":
    main()
