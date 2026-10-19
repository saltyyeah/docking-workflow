"""Rank experimental PDB candidates using transparent, limited heuristics.

This is a triage aid. Human origin, ligand pocket suitability, wild type,
assembly, completeness and non-standard residues must still be manually checked.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def number(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def score(row: dict) -> tuple:
    method = (row.get("experimental_method") or "").lower()
    resolution = number(row.get("resolution"), 99.0)
    coverage = number(row.get("coverage"), 0.0)
    return (
        "x-ray" in method,
        resolution <= 2.5,
        "diffraction" in method,
        coverage,
        -resolution,
        -number(row.get("pdbe_rank"), 9999.0),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidates", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--top", type=int, default=3)
    args = parser.parse_args()
    with args.candidates.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    groups: dict[str, list[dict]] = {}
    for row in rows:
        groups.setdefault(row["target"], []).append(row)
    selected = []
    for target, candidates in groups.items():
        for selection_rank, row in enumerate(sorted(candidates, key=score, reverse=True)[:args.top], start=1):
            row = dict(row)
            row["selection_rank"] = selection_rank
            row["manual_review"] = "true"
            row["review_notes"] = "Check species, assembly, mutations, binding ligand, pocket residues, cofactors and protonation."
            selected.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) + ["selection_rank", "manual_review", "review_notes"] if rows else []
    with args.output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(selected)
    print(f"selected={len(selected)} targets={len(groups)}")


if __name__ == "__main__":
    main()
