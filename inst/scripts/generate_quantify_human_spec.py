#!/usr/bin/env python3
"""Generate inst/extdata/quantify_human_immune_spec.tsv for human CellTypist labels.

Produces quantification rules for celltypist.Immune_All_High.cellclass and
celltypist.Immune_All_Low.cellclass metadata fields (CASCADE / human samples).

Usage (from package root):
  python3 inst/scripts/generate_quantify_human_spec.py \\
    --out inst/extdata/quantify_human_immune_spec.tsv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

GROUPING = "cDNA_ID"
DELIM = "__"
HIGH_FIELD = "celltypist.Immune_All_High.cellclass"
LOW_FIELD = "celltypist.Immune_All_Low.cellclass"

HIGH_CLASSES = [
    "B cells",
    "T cells",
    "Monocytes",
    "Macrophages",
    "DC",
    "ILC",
    "Plasma cells",
    "pDC",
    "Granulocytes",
    "Mast cells",
    "Erythrocytes",
    "Epithelial cells",
    "Endothelial cells",
    "Fibroblasts",
    "Cycling cells",
    "Unassigned",
    "Unknown",
    "Ambiguous",
]

LOW_CLASSES = [
    "Alveolar macrophages",
    "Classical monocytes",
    "Non-classical monocytes",
    "Macrophages",
    "DC1",
    "DC2",
    "pDC",
    "CD16+ NK cells",
    "NK cells",
    "Tcm/Naive helper T cells",
    "Tcm/Naive cytotoxic T cells",
    "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells",
    "Regulatory T cells",
    "MAIT cells",
    "gamma-delta T cells",
    "Memory B cells",
    "Naive B cells",
    "Plasma cells",
    "Unassigned",
    "Ambiguous",
    "Unknown",
]

LOW_HIGH_SCOPE = {
    "Alveolar macrophages": "Macrophages",
    "Classical monocytes": "Monocytes",
    "Non-classical monocytes": "Monocytes",
    "Macrophages": "Macrophages",
    "DC1": "DC",
    "DC2": "DC",
    "pDC": "pDC",
    "CD16+ NK cells": "ILC",
    "NK cells": "ILC",
    "Tcm/Naive helper T cells": "T cells",
    "Tcm/Naive cytotoxic T cells": "T cells",
    "Tem/Temra cytotoxic T cells": "T cells",
    "Tem/Trm cytotoxic T cells": "T cells",
    "Regulatory T cells": "T cells",
    "MAIT cells": "T cells",
    "gamma-delta T cells": "T cells",
    "Memory B cells": "B cells",
    "Naive B cells": "B cells",
    "Plasma cells": "Plasma cells",
}

LOW_T_CELL_LABELS = [
    "Tcm/Naive helper T cells",
    "Tcm/Naive cytotoxic T cells",
    "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells",
    "Regulatory T cells",
    "MAIT cells",
    "gamma-delta T cells",
]

LOW_NK_LABELS = ["CD16+ NK cells", "NK cells"]

LOW_SCORE_MODULES = [
    ("CytotoxicityScore", "Cytotoxicity_UCell"),
    ("InterferonResponseScore", "Interferon_Response_UCell"),
]

LOW_PCT_POSITIVE_GENES = ["PDCD1", "KLRK1", "GZMB"]  # subset of Quantify.PctPositive in R/Quantify.R

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


def build_spec_rows() -> list[str]:
    rows = [SPEC_HEADER]

    for cls in HIGH_CLASSES:
        rows.append(
            spec_row(join_target("ImmuneHigh", sanitize(cls)), HIGH_FIELD, cls, "sum")
        )

    for cls in LOW_CLASSES:
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

    for cls in LOW_T_CELL_LABELS:
        for score_suffix, score_field in LOW_SCORE_MODULES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cls), score_suffix),
                    LOW_FIELD,
                    cls,
                    "score",
                    score_field,
                    "quantiles",
                    HIGH_FIELD,
                    "T cells",
                )
            )

    for cls in LOW_NK_LABELS:
        for score_suffix, score_field in LOW_SCORE_MODULES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cls), score_suffix),
                    LOW_FIELD,
                    cls,
                    "score",
                    score_field,
                    "quantiles",
                    HIGH_FIELD,
                    "ILC",
                )
            )

    for cls in LOW_T_CELL_LABELS:
        for gene in LOW_PCT_POSITIVE_GENES:
            rows.append(
                spec_row(
                    join_target("ImmuneLow", sanitize(cls), "PctPositive", gene),
                    LOW_FIELD,
                    cls,
                    "pct_positive",
                    gene,
                    HIGH_FIELD,
                    "T cells",
                )
            )

    for cls in DIVERSITY_LOW_CLASSES:
        rows.append(
            spec_row(
                join_target("ImmuneLow", sanitize(cls), "Diversity"),
                LOW_FIELD,
                cls,
                "diversity",
                HIGH_FIELD,
                "T cells",
            )
        )

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("inst/extdata/quantify_human_immune_spec.tsv"),
    )
    args = parser.parse_args()

    rows = build_spec_rows()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"High classes ({len(HIGH_CLASSES)}): {HIGH_CLASSES}")
    print(f"Low classes ({len(LOW_CLASSES)}): {LOW_CLASSES}")
    print(f"Total spec rows (excl. header): {len(rows) - 1}")


if __name__ == "__main__":
    main()
