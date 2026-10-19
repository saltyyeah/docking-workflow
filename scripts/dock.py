"""Run Meeko preparation, native-ligand redocking, and query docking from CSV.

Required CSV columns: task_id,target,uniprot_id,pdb_id,receptor_prepared_pdb,
native_ligand_reference_sdf,native_ligand_prepared_sdf,query_ligand_prepared_sdf.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


REQUIRED = {"task_id", "target", "uniprot_id", "pdb_id", "receptor_prepared_pdb", "native_ligand_reference_sdf", "native_ligand_prepared_sdf", "query_ligand_prepared_sdf"}


def resolve(value: str, base: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def box_from_sdf(filename: Path, padding: float) -> dict:
    from rdkit import Chem
    molecule = Chem.SDMolSupplier(str(filename), removeHs=False)[0]
    if molecule is None or not molecule.GetNumConformers():
        raise RuntimeError(f"Invalid 3D SDF: {filename}")
    conformer = molecule.GetConformer()
    points = [conformer.GetAtomPosition(atom.GetIdx()) for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    if not points:
        raise RuntimeError(f"No ligand heavy atoms: {filename}")
    lower = [min(getattr(point, axis) for point in points) for axis in "xyz"]
    upper = [max(getattr(point, axis) for point in points) for axis in "xyz"]
    return {"center": [(a + b) / 2 for a, b in zip(lower, upper)], "size": [b - a + 2 * padding for a, b in zip(lower, upper)], "padding_A": padding}


def execute(command: list[str], workdir: Path, log, enabled: bool) -> None:
    log.write("$ " + " ".join(command) + "\n")
    if not enabled:
        return
    completed = subprocess.run(command, cwd=workdir, text=True, capture_output=True)
    log.write(completed.stdout); log.write(completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"Command failed: {command[0]}")


def process(row: dict, base: Path, output: Path, args) -> dict:
    summary = {key: row.get(key, "") for key in ("task_id", "target", "uniprot_id", "pdb_id")}
    summary["status"] = "failed"
    try:
        if REQUIRED - set(row):
            raise RuntimeError("Manifest has missing required columns")
        if row.get("manual_review", "").lower() in {"true", "1", "yes"}:
            raise RuntimeError("Row requires manual review")
        inputs = {name: resolve(row[name], base) for name in REQUIRED - {"task_id", "target", "uniprot_id", "pdb_id"}}
        missing = [name for name, value in inputs.items() if not value.is_file()]
        if missing:
            raise RuntimeError("Missing input: " + ", ".join(missing))
        task = output / row["task_id"]
        task.mkdir(parents=True, exist_ok=True)
        box = box_from_sdf(inputs["native_ligand_reference_sdf"], args.padding)
        (task / "box.json").write_text(json.dumps(box, indent=2), encoding="utf-8")
        config = task / "vina_box.txt"
        config.write_text("\n".join([f"receptor = receptor.pdbqt", *[f"center_{axis} = {value:.3f}" for axis, value in zip("xyz", box["center"])], *[f"size_{axis} = {value:.3f}" for axis, value in zip("xyz", box["size"])]] ) + "\n", encoding="utf-8")
        receptor, native, query = task / "receptor.pdbqt", task / "native.pdbqt", task / "query.pdbqt"
        preparation = [
            [args.meeko_receptor, "--read_pdb", str(inputs["receptor_prepared_pdb"]), "-o", str(task / "receptor"), "-p", "-v", "--box_center", *map(str, box["center"]), "--box_size", *map(str, box["size"])],
            [args.meeko_ligand, "-i", str(inputs["native_ligand_prepared_sdf"]), "-o", str(native)],
            [args.meeko_ligand, "-i", str(inputs["query_ligand_prepared_sdf"]), "-o", str(query)],
        ]
        vina = [args.vina_bin, "--receptor", str(receptor), "--config", str(config), "--scoring", "vina", "--exhaustiveness", str(args.exhaustiveness), "--num_modes", str(args.num_modes), "--energy_range", str(args.energy_range), "--seed", str(args.seed)]
        commands = preparation + [vina + ["--ligand", str(native), "--out", str(task / "native_redocked.pdbqt")], vina + ["--ligand", str(query), "--out", str(task / "query_docked.pdbqt")]]
        with (task / "commands.log").open("w", encoding="utf-8") as log:
            for command in commands:
                execute(command, task, log, args.execute)
        (task / "run_metadata.json").write_text(json.dumps({"input": row, "box": box, "commands": commands, "parameters": vars(args)}, indent=2), encoding="utf-8")
        summary.update({"status": "completed" if args.execute else "planned", "box_center": ";".join(f"{x:.3f}" for x in box["center"]), "box_size": ";".join(f"{x:.3f}" for x in box["size"]), "task_dir": str(task)})
    except Exception as error:
        summary["reason"] = str(error)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path); parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute", action="store_true"); parser.add_argument("--padding", type=float, default=5.0)
    parser.add_argument("--seed", type=int, default=2026); parser.add_argument("--exhaustiveness", type=int, default=32); parser.add_argument("--num-modes", type=int, default=9); parser.add_argument("--energy-range", type=float, default=3.0)
    parser.add_argument("--vina-bin", default="vina"); parser.add_argument("--meeko-receptor", default="mk_prepare_receptor.py"); parser.add_argument("--meeko-ligand", default="mk_prepare_ligand.py")
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle: rows = list(csv.DictReader(handle))
    results = [process(row, args.manifest.resolve().parent, args.output, args) for row in rows]
    fields = sorted({key for row in results for key in row})
    for filename, report in (("results.csv", results), ("failed_tasks.csv", [row for row in results if row["status"] == "failed"])):
        with (args.output / filename).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(report)
    if any(row["status"] == "failed" for row in results): sys.exit(1)


if __name__ == "__main__":
    main()
