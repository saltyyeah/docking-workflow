"""Extract one crystallographic ligand instance as SDF via RCSB ModelServer.

The input CSV must identify an individual ligand instance with pdb_id,
auth_asym_id,auth_seq_id,comp_id. The output SDF preserves experimental
coordinates and must be retained unchanged as the reference ligand for Box
calculation and RMSD evaluation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


REQUIRED = {"pdb_id", "auth_asym_id", "auth_seq_id", "comp_id"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ligands", type=Path, help="CSV of individual ligand instances")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.ligands.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not REQUIRED.issubset(rows[0]):
        raise SystemExit("Input needs pdb_id,auth_asym_id,auth_seq_id,comp_id columns.")
    args.output.mkdir(parents=True, exist_ok=True)
    report = []
    for row in rows:
        pdb_id = row["pdb_id"].strip().lower()
        params = urlencode({"auth_asym_id": row["auth_asym_id"].strip(), "auth_seq_id": row["auth_seq_id"].strip(), "encoding": "sdf"})
        filename = f"{pdb_id}_{row['auth_asym_id'].strip()}_{row['auth_seq_id'].strip()}_{row['comp_id'].strip()}.sdf"
        destination = args.output / filename
        try:
            url = f"https://models.rcsb.org/v1/{pdb_id}/ligand?{params}"
            request = Request(url, headers={"User-Agent": "vina-project/1.0"})
            with urlopen(request, timeout=60) as response:
                destination.write_bytes(response.read())
            report.append({**row, "reference_sdf": str(destination), "status": "extracted", "reason": ""})
        except (HTTPError, URLError, TimeoutError) as error:
            report.append({**row, "reference_sdf": "", "status": "failed", "reason": str(error)})
    fields = list(rows[0]) + ["reference_sdf", "status", "reason"]
    with (args.output / "ligand_extract_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(report)
    print(f"extracted={sum(row['status'] == 'extracted' for row in report)}")


if __name__ == "__main__":
    main()
