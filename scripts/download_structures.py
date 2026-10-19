"""Download selected PDB coordinate files from RCSB by PDB ID."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selected", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.selected.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for row in rows:
        pdb_id = row.get("pdb_id", "").strip().lower()
        destination = args.output / f"{pdb_id}.pdb"
        try:
            request = Request(f"https://files.rcsb.org/download/{pdb_id}.pdb", headers={"User-Agent": "vina-project/1.0"})
            with urlopen(request, timeout=60) as response:
                destination.write_bytes(response.read())
            results.append({"pdb_id": pdb_id, "file": str(destination), "status": "downloaded", "reason": ""})
        except (HTTPError, URLError, TimeoutError) as error:
            results.append({"pdb_id": pdb_id, "file": "", "status": "failed", "reason": str(error)})
    with (args.output / "download_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("pdb_id", "file", "status", "reason")); writer.writeheader(); writer.writerows(results)
    print(f"downloaded={sum(row['status'] == 'downloaded' for row in results)}")


if __name__ == "__main__":
    main()
