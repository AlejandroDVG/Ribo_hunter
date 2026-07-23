#!/bin/bash
#SBATCH --job-name=ribohunter_dl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# OPTIONAL: SLURM array wrapper for download.sh -- only needed if you're running
# this on a Slurm cluster. For a single machine, just run download.sh directly.
#
# Submit with, e.g.:
#   mkdir -p logs
#   sbatch --array=1-N%5 download.slurm manifest.csv /path/to/outdir
#   (%5 = 5 concurrent downloads at a time -- start conservative on a shared cluster)
#
# Before submitting the full array, run one task end to end first and check the
# result: correct byte count, `gzip -t` passes, real FASTQ content on the first
# line. Confirming one sample costs a minute; discovering a path bug 500 tasks
# into a real run does not.

set -euo pipefail

MANIFEST="${1:?Usage: sbatch --array=1-N download.slurm MANIFEST.csv OUTDIR}"
OUTDIR="${2:?Usage: sbatch --array=1-N download.slurm MANIFEST.csv OUTDIR}"

# Slurm copies this script into a per-job spool directory before running it, so a
# BASH_SOURCE-relative lookup of the sibling download.sh would resolve to the
# spool copy's location, not where you submitted from -- use SLURM_SUBMIT_DIR.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR not set -- run via sbatch, not directly}"
"$SCRIPT_DIR/download.sh" "$MANIFEST" "$OUTDIR" --array-mode
