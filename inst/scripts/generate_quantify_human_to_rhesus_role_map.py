#!/usr/bin/env python3
"""Generate inst/extdata/quantify_human_to_rhesus_role_map.tsv.

Maps CellTypist Immune_All High/Low labels to RIRA metadata fields for
coerceToRIRA quantification. Prefer Low labels; High is fallback when Low
is absent or unmapped.

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

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generate_quantify_human_spec import (  # noqa: E402
    HARDCODED_HIGH,
    HARDCODED_LOW,
    LOW_HIGH_SCOPE,
)

IMMUNE_FIELD = "RIRA_Immune_v2.cellclass"
T_AND_NK_FIELD = "RIRA_TNK_v2.cellclass"
MYELOID_COARSE_FIELD = "RIRA_Myeloid_v3.coarseclass"
MYELOID_FINE_FIELD = "RIRA_Myeloid_v3.cellclass"

HIGH_SOURCE = "Immune_All_High"
LOW_SOURCE = "Immune_All_Low"

# High-level Immune coarse mapping (CellTypist High onto RIRA Immune).
# Cycling cells High stays T_NK/T as a coarse fallback when Low is absent.
HIGH_IMMUNE_MAP = {
    "B cells": ("Bcell", "B"),
    "B-cell lineage": ("Bcell", "B"),
    "Plasma cells": ("Bcell", "B"),
    "T cells": ("T_NK", "T"),
    "Cycling cells": ("T_NK", "T"),
    "ILC": ("T_NK", "NK_ILC"),
    "ILC precursor": ("T_NK", "NK_ILC"),
    "Double-negative thymocytes": ("T_NK", "T"),
    "Double-positive thymocytes": ("T_NK", "T"),
    "ETP": ("T_NK", "T"),
    "Macrophages": ("Myeloid", "Myeloid_APC"),
    "Monocytes": ("Myeloid", "Myeloid_APC"),
    "DC": ("Myeloid", "Myeloid_APC"),
    "DC precursor": ("Myeloid", "Myeloid_APC"),
    "MNP": ("Myeloid", "Myeloid_APC"),
    "Mono-mac": ("Myeloid", "Myeloid_APC"),
    "Monocyte precursor": ("Myeloid", "Myeloid_APC"),
    "Myelocytes": ("Myeloid", "Myeloid_APC"),
    "Promyelocytes": ("Myeloid", "Myeloid_APC"),
    "Megakaryocytes/platelets": ("Myeloid", "Myeloid_APC"),
    "Megakaryocyte precursor": ("Myeloid", "Myeloid_APC"),
    "Early MK": ("Myeloid", "Myeloid_APC"),
    "Mast cells": ("Myeloid", "Myeloid_APC"),
    "Granulocytes": ("Myeloid", "Myeloid_APC"),
    "pDC": ("Myeloid", "Myeloid_APC"),
    "pDC precursor": ("Myeloid", "Myeloid_APC"),
    "Epithelial cells": ("Epithelial", "Epithelial"),
    "Erythrocytes": ("Erythrocyte", "Erythrocyte"),
    "Erythroid": ("Erythrocyte", "Erythrocyte"),
    "Endothelial cells": ("Stromal", "Stromal"),
    "Fibroblasts": ("Stromal", "Stromal"),
    "HSC/MPP": ("Non-Immune", "HSC"),
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

# NKT are T cells (CellTypist High=T cells); CD8-associated markers → CD8_like.
CD8_LIKE = {
    "Tcm/Naive cytotoxic T cells",
    "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells",
    "Trm cytotoxic T cells",
    "CD8a/a",
    "CD8a/b(entry)",
    "MAIT cells",
    "NKT cells",
}

GAMMA_DELTA_LIKE = {
    "gamma-delta T cells",
    "CRTAM+ gamma-delta T cells",
    "Cycling gamma-delta T cells",
}

# True NK Low labels → RIRA TNK "NK Cells".
NK_LIKE = {
    "NK cells",
    "CD16+ NK cells",
    "CD16- NK cells",
    "Cycling NK cells",
    "Transitional NK",
}

# Helper ILCs (non-cytotoxic / Th-analog) → RIRA TNK "Other" (no ILC class).
# Role remains NK_ILC; Immune coarse remains T_NK.
HELPER_ILC_LIKE = {
    "ILC",
    "ILC precursor",
    "ILC1",
    "ILC2",
    "ILC3",
}

# Cycling Low labels: map by subtype (do not inherit High Cycling cells → T).
# Values: (immune_label, role, optional_tnk_label, optional_myeloid_coarse)
LOW_CYCLING_MAP = {
    "Cycling B cells": ("Bcell", "B", None, None),
    "Cycling T cells": ("T_NK", "T", "Other", None),
    "Cycling gamma-delta T cells": ("T_NK", "GD", "Gamma/Delta Cells", None),
    "Cycling NK cells": ("T_NK", "NK_ILC", "NK Cells", None),
    "Cycling DCs": ("Myeloid", "Myeloid_APC", None, "DC"),
    "Cycling monocytes": ("Myeloid", "Myeloid_APC", None, "Monocytes"),
}

LOW_T_AND_NK_MAP = {}
for label in CD4_LIKE:
    LOW_T_AND_NK_MAP[label] = ("CD4+ T Cells", "CD4_like")
for label in CD8_LIKE:
    LOW_T_AND_NK_MAP[label] = ("CD8+ T Cells", "CD8_like")
for label in GAMMA_DELTA_LIKE:
    LOW_T_AND_NK_MAP[label] = ("Gamma/Delta Cells", "GD")
for label in NK_LIKE:
    LOW_T_AND_NK_MAP[label] = ("NK Cells", "NK_ILC")
for label in HELPER_ILC_LIKE:
    LOW_T_AND_NK_MAP[label] = ("Other", "NK_ILC")

LOW_MYELOID_COARSE = {
    "Classical monocytes": "Monocytes",
    "Non-classical monocytes": "Monocytes",
    "Cycling monocytes": "Monocytes",
    "Monocytes": "Monocytes",
    "DC": "DC",
    "DC1": "DC",
    "DC2": "DC",
    "DC3": "DC",
    "Migratory DCs": "DC",
    "Transitional DC": "DC",
    "Cycling DCs": "DC",
    "Macrophages": "Macrophages",
    "Alveolar macrophages": "Alv. mac.",
    "Erythrophagocytic macrophages": "Macrophages",
    "Intermediate macrophages": "Macrophages",
    "Intestinal macrophages": "Macrophages",
    "Kidney-resident macrophages": "Macrophages",
    "Kupffer cells": "Macrophages",
    "Hofbauer cells": "Macrophages",
    "pDC": "pDC",
    "Myelocytes": "Myelocytes",
    "Promyelocytes": "Promyelocytes",
    "Neutrophils": "Macrophages",
}

LOW_MYELOID_FINE = {
    "Classical monocytes": "CD14+ Monocytes",
    "Non-classical monocytes": "CD16+ Monocytes",
    "Inflammatory Monocytes": "Inflammatory Monocytes",
}

ILC_PARENTS = {"ILC", "ILC precursor"}


def immune_from_low_parent(low_label: str) -> tuple[str, str] | None:
    parent = LOW_HIGH_SCOPE.get(low_label)
    if not parent:
        return None
    return HIGH_IMMUNE_MAP.get(parent)


def build_rows() -> list[str]:
    header = "\t".join(
        [
            "humanSourceField",
            "humanLabel",
            "rhesusTargetField",
            "rhesusLabel",
            "lineageRole",
        ]
    )
    rows = [header]
    seen: set[tuple[str, str, str, str]] = set()

    def add(source: str, label: str, target: str, rhesus: str, role: str) -> None:
        key = (source, label, target, rhesus)
        if key in seen:
            return
        seen.add(key)
        rows.append("\t".join([source, label, target, rhesus, role]))

    for high in HARDCODED_HIGH:
        mapped = HIGH_IMMUNE_MAP.get(high)
        if mapped:
            add(HIGH_SOURCE, high, IMMUNE_FIELD, mapped[0], mapped[1])

    for low in HARDCODED_LOW:
        parent = LOW_HIGH_SCOPE.get(low)

        if low in LOW_CYCLING_MAP:
            immune_label, role, tnk_label, myeloid_coarse = LOW_CYCLING_MAP[low]
            add(LOW_SOURCE, low, IMMUNE_FIELD, immune_label, role)
            if tnk_label:
                add(LOW_SOURCE, low, T_AND_NK_FIELD, tnk_label, role)
            if myeloid_coarse:
                add(LOW_SOURCE, low, MYELOID_COARSE_FIELD, myeloid_coarse, role)
            continue

        # T cells, helper ILCs, and NK Lows parented under CellTypist High ILC.
        if (
            parent in ({"T cells"} | ILC_PARENTS) or low in HELPER_ILC_LIKE
        ) and low in LOW_T_AND_NK_MAP:
            t_and_nk_label, role = LOW_T_AND_NK_MAP[low]
            add(LOW_SOURCE, low, IMMUNE_FIELD, "T_NK", role)
            add(LOW_SOURCE, low, T_AND_NK_FIELD, t_and_nk_label, role)
        elif parent in {"B cells", "B-cell lineage", "Plasma cells"}:
            add(LOW_SOURCE, low, IMMUNE_FIELD, "Bcell", "B")
        elif parent in {
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
        }:
            immune = immune_from_low_parent(low) or ("Myeloid", "Myeloid_APC")
            add(LOW_SOURCE, low, IMMUNE_FIELD, immune[0], immune[1])
            coarse = LOW_MYELOID_COARSE.get(low)
            if coarse:
                add(LOW_SOURCE, low, MYELOID_COARSE_FIELD, coarse, "Myeloid_APC")
            fine = LOW_MYELOID_FINE.get(low)
            if fine:
                add(LOW_SOURCE, low, MYELOID_FINE_FIELD, fine, "Myeloid_APC")
        else:
            immune = immune_from_low_parent(low)
            if immune:
                add(LOW_SOURCE, low, IMMUNE_FIELD, immune[0], immune[1])

    return rows


def main() -> None:
    out = Path("inst/extdata/quantify_human_to_rhesus_role_map.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows) - 1} map rows to {out}")


if __name__ == "__main__":
    main()
