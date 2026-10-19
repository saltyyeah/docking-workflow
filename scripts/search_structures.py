"""Search PDBe best experimental structures for UniProt target IDs.

Input CSV must contain: target,uniprot_id. Output is a long candidate table;
it does not claim that any hit is a suitable docking receptor.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

API = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/"
FIELDS = ("target", "uniprot_id", "pdbe_rank", "pdb_id", "chain_id", "coverage", "resolution", "experimental_method", "unp_start", "unp_end")


def fetch(target: str, uniprot_id: str, timeout: int) -> tuple[list[dict], dict | None]:
    uniprot_id = uniprot_id.strip().upper()
    try:
        request = Request(API + uniprot_id, headers={"User-Agent": "vina-project/1.0"})
        with urlopen(request, timeout=timeout) as response:
            data = json.load(response)
        structures = data.get(uniprot_id, [])
        if not structures:
            return [], {"target": target, "uniprot_id": uniprot_id, "reason": "No PDBe mapping"}
        rows = []
        for rank, item in enumerate(structures, start=1):
            rows.append({"target": target, "uniprot_id": uniprot_id, "pdbe_rank": rank, "pdb_id": item.get("pdb_id"), "chain_id": item.get("chain_id"), "coverage": item.get("coverage"), "resolution": item.get("resolution"), "experimental_method": item.get("experimental_method"), "unp_start": item.get("unp_start", item.get("start")), "unp_end": item.get("unp_end", item.get("end"))})
        return rows, None
    except (HTTPError, URLError, TimeoutError, ValueError, json.JSONDecodeError) as error:
        return [], {"target": target, "uniprot_id": uniprot_id, "reason": str(error)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()
    with args.targets.open(encoding="utf-8-sig", newline="") as handle:
        targets = list(csv.DictReader(handle))
    if not targets or not {"target", "uniprot_id"}.issubset(targets[0]):
        raise SystemExit("Input needs target and uniprot_id columns.")
    rows, failures = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(fetch, row["target"], row["uniprot_id"], args.timeout) for row in targets]
        for future in as_completed(futures):
            found, failed = future.result()
            rows.extend(found)
            if failed:
                failures.append(failed)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "pdb_candidates.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS); writer.writeheader(); writer.writerows(sorted(rows, key=lambda x: (x["target"], x["pdbe_rank"])))
    with (args.output / "pdb_search_failed.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("target", "uniprot_id", "reason")); writer.writeheader(); writer.writerows(failures)
    print(f"candidates={len(rows)} failures={len(failures)}")


if __name__ == "__main__":
    main()
