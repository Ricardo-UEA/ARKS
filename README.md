# ARKS

Ancestry-Representative K-mer Set — a contamination-screened, multi-ancestry pan-genome k-mer database for host-read depletion in clinical metagenomics. This repository contains the pipeline and analysis code used to build and benchmark ARKS.

## Overview

ARKS integrates human reference datasets — GRCh37, GRCh38, T2T-CHM13+Y, pan-genome assemblies for African, Arab, Chinese, and Khoe-San populations, an additional pan-genome reference, and IPD-IMGT/HLA — into a single k-mer database, screened for microbial contamination and benchmarked against single-reference depletion (GRCh38/T2T) using both CAMISIM synthetic metagenomes and direct read-depletion comparisons.

This repository holds **code only**. Raw sequence data, intermediate assemblies, and pipeline outputs are not included — see [Data & outputs](#data--outputs).

## Repository structure

### `references/`
K-mer construction for the baseline references (GRCh37, GRCh38, T2T-CHM13+Y): Jellyfish counting (`jellyfish_references.sh`), Kraken classification, Snakemake driver (`00_snakemake.sh`).

### `african_pan/`, `arab_pan_genome/`, `chinese_pan_genome/`, `khoisan_pan_genome/`, `pan_genome_2/`
Per-dataset pan-genome construction pipelines, each following the same pattern: BCALM2 assembly (`bcalm_assembly.sh` / `bcalm2_kmer_assembly.sh`) → FASTA/FASTQ conversion → post-assembly validation → Kraken classification, driven by a numbered Snakemake wrapper (`0N_snakemake.sh`). Some folders also contain convenience variants used during development (`full_script_simple.sh`, `quick_classify.sh`, `kraken_quick.sh`).

### `ipd_imgt_hla_db/`
Same construction pattern applied to the IPD-IMGT/HLA reference.

### `unique_kmers/`
Contamination screening, split into one subfolder per dataset (`GRCh37/`, `GRCh38/`, `T2T/`, each pan-genome, and `ipd_hla_immune/`): FASTA conversion, Kraken classification, and microbial k-mer removal (`microbial_removal.py`).

### `all_kmers/`
Merges all per-dataset k-mer sets into the final ARKS database: Jellyfish merge (`all_kmers_jellyfish.sh`), DuckDB conversion, k-mer summary table generation, UpSet plot generation (`upset_plot_final.py`), and validation/stats.

### `african_sydney/`
Host-depletion benchmarking: BBDuk depletion runs comparing ARKS, T2T, and the legacy PanHuman database, each followed by Kraken classification of the depleted reads (before/after).

### `simulations/`
CAMISIM-based synthetic metagenome benchmarking.
- **`simulations/CAMISIM/`** — vendored copy of [CAMISIM](https://github.com/CAMI-challenge/CAMISIM) (legacy Python standalone, pre-2.0/Nextflow). Third-party code — see [Third-party code](#third-party-code).
- **`simulations/CAMISIM/meta_simulations_2026_v2/`** — our own simulation run built on top of CAMISIM: config generation (`generate_camisim_configs.py`), BBDuk depletion, Kraken classification, read-count/assembly evaluation, and ground-truth summary generation. This is the CAMISIM benchmarking pipeline reported in the paper.

## Third-party code

`simulations/CAMISIM/` (excluding `meta_simulations_2026_v2/`) is vendored from [CAMI-challenge/CAMISIM](https://github.com/CAMI-challenge/CAMISIM), licensed under Apache-2.0, and is not original to this project. If reusing it, cite:

> Fritz\*, Hofmann\*, et al. (2019). CAMISIM: Simulating metagenomes and microbial communities. *Microbiome*, 7:17. doi:[10.1186/s40168-019-0633-6](https://doi.org/10.1186/s40168-019-0633-6)

## Dependencies

Snakemake, BCALM2, Jellyfish, KrakenUniq/Kraken2, BBTools (BBDuk), CAMISIM, DuckDB, Singularity, SLURM (developed on the Ada HPC cluster).

## Data & outputs

This repository tracks scripts only. Within each folder, subdirectories such as `output_files/`, `filtered_assembly/`, `input_files/`, `kraken_reports/`, `stats_files/`, `job_errors/`, and `.snakemake/` hold generated data, intermediate assemblies, and logs on the working copy — only the script files from these subdirectories are included here. The ARKS k-mer database and raw sequencing data are available at [Zenodo DOI — TBD] / on request.

## Citation

If you use this pipeline or the ARKS database, please cite:

> [Citation — TBD pending publication]

## License

MIT (see `LICENSE`), except for vendored code under `simulations/CAMISIM/` (Apache-2.0 — see [Third-party code](#third-party-code)).
