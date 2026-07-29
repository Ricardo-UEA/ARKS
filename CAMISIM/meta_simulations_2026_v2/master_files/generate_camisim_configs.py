#!/usr/bin/env python3
"""
Generate CAMISIM configs + per-sample id_to_genome and metadata subsets.

Expected input files:
  - camisim_config_template.ini (template config file)
  - id_to_genome_master.tsv     (2 columns, no header: genome_ID <tab> genome_path)
  - metadata_master.tsv         (TSV with header: genome_ID, OTU, NCBI_ID, novelty_category)

Outputs created:
  - output_files/meta_sample_X/           (CAMISIM output_directory target)
  - meta_simulations_2026/
      - configs/config_X.ini
      - id_to_genomes/id_to_genome_X.tsv
      - metadata/metadata_X.tsv
"""

from __future__ import annotations

import random
import re
from pathlib import Path
from typing import FrozenSet, List, Set

import pandas as pd


# ============================================================
# ======================= USER CONFIG ========================
# ============================================================

# Base directory on the HPC
BASE_DIR = Path("/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2")

# Inputs
TEMPLATE_INI = BASE_DIR / "master_files/camisim_config_master.ini"
ID_TO_GENOME_MASTER = BASE_DIR / "master_files/id_to_genome_master.tsv"
METADATA_MASTER = BASE_DIR / "master_files/metadata_master.tsv"

# How many configs/samples to generate
N_SAMPLES = 100

# Randomly choose K genomes per sample, where K is between MIN_GENOMES and MAX_GENOMES
MIN_GENOMES = 55
MAX_GENOMES = 65

# Reproducibility
RANDOM_SEED = 22

# If True, avoid identical genome sets across samples
UNIQUE_COMBOS = False

# Output layout
OUTPUT_ROOT = BASE_DIR / "output_files"                    # CAMISIM output: output_files/meta_sample_X
META_SIM_ROOT = BASE_DIR / "meta_simulations_2026"         # Supporting files root
CONFIG_OUTDIR = BASE_DIR / "configs"                  # configs/config_X.ini
ID_TO_GENOME_OUTDIR = BASE_DIR / "id_to_genomes"      # id_to_genomes/id_to_genome_X.tsv
METADATA_OUTDIR = BASE_DIR / "metadata"               # metadata/metadata_X.tsv

# ============================================================
# ============================================================


def ensure_dir(p: Path) -> None:
    """Create directory if it doesn't exist."""
    p.mkdir(parents=True, exist_ok=True)


def read_id_to_genome(path: Path) -> pd.DataFrame:
    """Read id_to_genome file (2 columns, no header)."""
    df = pd.read_csv(path, sep="\t", header=None, names=["genome_ID", "genome_path"], dtype=str)
    if df.empty:
        raise ValueError(f"id_to_genome_master is empty: {path}")
    if df["genome_ID"].isna().any():
        raise ValueError(f"Found missing genome_ID values in: {path}")
    return df


def read_metadata(path: Path) -> pd.DataFrame:
    """Read metadata file (TSV with header)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        raise ValueError(f"metadata_master is empty: {path}")
    if "genome_ID" not in df.columns:
        raise ValueError(f"metadata_master must contain a 'genome_ID' column. Found: {list(df.columns)}")
    return df


def patch_ini_template(
    template_text: str,
    output_directory: Path,
    metadata_path: Path,
    id_to_genome_path: Path,
    genomes_total: int,
) -> str:
    """
    Patch the config template with sample-specific values.
    Updates in [Main]: output_directory
    Updates in [community0]: metadata, id_to_genome_file, genomes_total, num_real_genomes
    """
    replacements = {
        r"^output_directory=.*$": f"output_directory={output_directory.as_posix()}",
        r"^metadata=.*$": f"metadata={metadata_path.as_posix()}",
        r"^id_to_genome_file=.*$": f"id_to_genome_file={id_to_genome_path.as_posix()}",
        r"^genomes_total=.*$": f"genomes_total={genomes_total}",
        r"^num_real_genomes=.*$": f"num_real_genomes={genomes_total}",
    }

    out_lines = []
    for line in template_text.splitlines(True):
        stripped = line.strip()
        replaced = False
        for pat, repl in replacements.items():
            if re.match(pat, stripped):
                out_lines.append(repl + ("\n" if line.endswith("\n") else ""))
                replaced = True
                break
        if not replaced:
            out_lines.append(line)

    return "".join(out_lines)


def renumber_genomes(id_sub: pd.DataFrame, meta_sub: pd.DataFrame) -> tuple:
    """
    Renumber genome IDs sequentially (Genome1, Genome2, ...) for the subset.
    Returns renumbered dataframes.
    """
    # Create mapping from old ID to new sequential ID
    old_ids = sorted(id_sub["genome_ID"].tolist(), key=lambda x: int(x.replace("Genome", "")))
    id_mapping = {old_id: f"Genome{i+1}" for i, old_id in enumerate(old_ids)}
    
    # Apply mapping
    id_sub = id_sub.copy()
    meta_sub = meta_sub.copy()
    
    id_sub["genome_ID"] = id_sub["genome_ID"].map(id_mapping)
    meta_sub["genome_ID"] = meta_sub["genome_ID"].map(id_mapping)
    
    # Update OTU to match new numbering
    meta_sub["OTU"] = meta_sub["genome_ID"].apply(lambda x: x.replace("Genome", ""))
    
    return id_sub, meta_sub


def main() -> None:
    # Check inputs exist
    for p in (TEMPLATE_INI, ID_TO_GENOME_MASTER, METADATA_MASTER):
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p}")

    if MIN_GENOMES > MAX_GENOMES:
        raise ValueError("MIN_GENOMES cannot be greater than MAX_GENOMES")

    random.seed(RANDOM_SEED)

    # Load inputs
    id_df = read_id_to_genome(ID_TO_GENOME_MASTER)
    meta_df = read_metadata(METADATA_MASTER)

    # Keep only metadata rows that have genome_ID present in id_to_genome_master
    pool_ids = set(id_df["genome_ID"].tolist())
    meta_df = meta_df[meta_df["genome_ID"].isin(pool_ids)].copy()

    if meta_df.empty:
        raise ValueError(
            "After filtering metadata_master to IDs present in id_to_genome_master, no rows remain. "
            "Check that genome_ID values match between the two files."
        )

    genome_ids: List[str] = sorted(pool_ids, key=lambda x: int(x.replace("Genome", "")))
    n_pool = len(genome_ids)

    if MAX_GENOMES > n_pool:
        raise ValueError(
            f"MAX_GENOMES={MAX_GENOMES} but your pool only has {n_pool} genomes. "
            f"Lower MAX_GENOMES or add more genomes to the master files."
        )

    # Prepare output dirs
    ensure_dir(OUTPUT_ROOT)
    ensure_dir(META_SIM_ROOT)
    ensure_dir(CONFIG_OUTDIR)
    ensure_dir(ID_TO_GENOME_OUTDIR)
    ensure_dir(METADATA_OUTDIR)

    template_text = TEMPLATE_INI.read_text()

    seen: Set[FrozenSet[str]] = set()

    # Generate samples
    print(f"Generating {N_SAMPLES} samples...")
    for i in range(1, N_SAMPLES + 1):
        k = random.randint(MIN_GENOMES, MAX_GENOMES)

        # Pick genome set (optionally unique)
        attempts = 0
        max_attempts = 10000
        while True:
            chosen = frozenset(random.sample(genome_ids, k))
            if not UNIQUE_COMBOS or chosen not in seen:
                seen.add(chosen)
                break
            attempts += 1
            if attempts > max_attempts:
                raise RuntimeError(
                    f"Could not find unique combination after {max_attempts} attempts. "
                    f"Consider reducing N_SAMPLES or setting UNIQUE_COMBOS=False."
                )

        # Subset dataframes
        id_sub = id_df[id_df["genome_ID"].isin(chosen)].copy()
        meta_sub = meta_df[meta_df["genome_ID"].isin(chosen)].copy()

        # Renumber genomes sequentially (Genome1, Genome2, ...)
        id_sub, meta_sub = renumber_genomes(id_sub, meta_sub)

        # Sort by genome number
        id_sub = id_sub.sort_values("genome_ID", key=lambda x: x.str.replace("Genome", "").astype(int))
        meta_sub = meta_sub.sort_values("genome_ID", key=lambda x: x.str.replace("Genome", "").astype(int))

        # Output paths
        id_out = ID_TO_GENOME_OUTDIR / f"id_to_genome_{i}.tsv"
        meta_out = METADATA_OUTDIR / f"metadata_{i}.tsv"
        cfg_out = CONFIG_OUTDIR / f"config_{i}.ini"
        out_dir = OUTPUT_ROOT / f"meta_sample_{i}"

        # Write id_to_genome (no header)
        id_sub.to_csv(id_out, sep="\t", index=False, header=False)
        
        # Write metadata (with header)
        meta_sub.to_csv(meta_out, sep="\t", index=False)

        # Patch config and write
        cfg_text = patch_ini_template(
            template_text=template_text,
            output_directory=out_dir,
            metadata_path=meta_out,
            id_to_genome_path=id_out,
            genomes_total=len(chosen),
        )
        cfg_out.write_text(cfg_text)

        if i % 10 == 0:
            print(f"  Generated {i}/{N_SAMPLES} samples...")

    print("\n✅ Done!")
    print(f"Pool size: {n_pool} genomes")
    print(f"Generated {N_SAMPLES} samples with K in [{MIN_GENOMES}, {MAX_GENOMES}]")
    print(f"\nOutput locations:")
    print(f"  Configs:        {CONFIG_OUTDIR}")
    print(f"  id_to_genomes:  {ID_TO_GENOME_OUTDIR}")
    print(f"  metadata:       {METADATA_OUTDIR}")
    print(f"  CAMISIM output: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
