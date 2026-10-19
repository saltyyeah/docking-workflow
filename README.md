# Protein structure to Vina project

A reusable Python workflow beginning with UniProt IDs: find experimental PDB structures, produce a reviewable shortlist, download reviewed structures, extract a chosen co-crystal ligand instance, and run Meeko plus AutoDock Vina.

## Install

```bash
conda env create -f environment.yml
conda activate wxy_vina_meeko_py311
```

Validated versions: Python 3.11, Meeko 0.8.0, AutoDock Vina 1.2.7, and RDKit 2024.09.6.

## Workflow

### 1. Search experimental structures

Provide a UTF-8 CSV containing `target,uniprot_id`.

```bash
python scripts/search_structures.py targets.csv --output data/search
```

The output includes PDBe rank, chain, UniProt coverage, experimental method, resolution, and UniProt residue span. It is a candidate list, not a docking decision.

### 2. Rank and manually review candidates

```bash
python scripts/rank_structures.py data/search/pdb_candidates.csv \
  --output data/selected_candidates.csv --top 3
```

The heuristic prefers X-ray structures, resolution <=2.5 A, coverage, then PDBe rank. All rows are intentionally set to `manual_review=true`. Before docking, review human origin, wild type, biological assembly, chains, bound drug-like ligand, pocket residues, missing residues, metals/cofactors, and protonation.

### 3. Download reviewed PDB structures

```bash
python scripts/download_structures.py data/selected_candidates.csv --output data/pdb
```

### 4. Extract experimental ligand coordinates

Create `ligand_instances.csv` with exactly one chosen ligand instance per row:

```text
pdb_id,auth_asym_id,auth_seq_id,comp_id
```

Then run:

```bash
python scripts/extract_ligand.py ligand_instances.csv --output data/ligands
```

The script calls RCSB ModelServer with `encoding=sdf` to obtain the specific ligand instance in its experimental coordinate frame. Keep this file unchanged as `native_ligand_reference_sdf`. Do not use ideal CCD coordinates as the reference ligand. Create a separately chemistry-reviewed SDF, retaining the reference heavy-atom coordinates, for `native_ligand_prepared_sdf`.

### 5. Create and run a docking manifest

Required columns:

```text
task_id,target,uniprot_id,pdb_id,receptor_prepared_pdb,native_ligand_reference_sdf,native_ligand_prepared_sdf,query_ligand_prepared_sdf
```

The receptor PDB must be repaired, protonated, and stripped of the chosen reference ligand before it enters the manifest. Run a dry plan first, then execute after review:

```bash
python scripts/dock.py docking_manifest.csv --output output/plan
python scripts/dock.py docking_manifest.csv --output output/run --execute
```

Each task calculates its own Box from reference ligand heavy atoms:

```text
center_i = (min_i + max_i) / 2
size_i = (max_i - min_i) + 2 * padding
```

Defaults are 5 A padding, Vina scoring, exhaustiveness 32, 9 poses, energy range 3 kcal/mol, and seed 2026. Keep parameters consistent within a comparative batch.

## Outputs and QC

Per task, the docking runner writes Box JSON/text, PDBQT files, docking poses, `commands.log`, and metadata. Batch files are `results.csv` and `failed_tasks.csv`.

Redock the native ligand with the same receptor and Box as the query ligand. Export the result using `mk_export.py`, calculate symmetry-aware heavy-atom RMSD against the unchanged reference SDF, and investigate failed redocking before screening. A common top-pose threshold is RMSD <=2.0 A.

## Example templates

`examples/1iep/targets.csv` is a minimal PDBe search input for ABL1. `examples/1iep/docking_manifest.csv` shows the manifest layout; no coordinate files are committed. Add only public or permitted example coordinates when testing.

## References

- [PDBe best structures mapping API](https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/)
- [RCSB ModelServer API](https://models.rcsb.org/)
- [AutoDock Vina basic docking](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)
- [Meeko tutorial](https://github.com/forlilab/Meeko/blob/develop/docs/source/tutorial1.rst)
