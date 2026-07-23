#!/bin/bash
# Download FASTQ files listed in a ribohunter manifest (build_manifest.py output).
#
# Manifest schema: one row per GSM, not per run (see build_manifest.py's docstring
# for why: one-row-per-run let two downloader tasks race on the same output filename
# for a multi-run GSM, silently losing data). A GSM with more than one SRA run lists
# all of them in `runs`/`fastq_ftp_by_run`, pipe-separated per run, semicolon-separated
# per file within a run.
#
# Method: direct FTP download from ENA (fastq_ftp_by_run column), not SRA-toolkit
#   prefetch/fasterq-dump -- ENA already serves gzipped FASTQ directly. For a
#   multi-run GSM, each run is downloaded to its own temp file and then
#   concatenated (gzip streams concatenate validly) into the single final file
#   per read position, in run order. Output layout:
#   <OUTDIR>/<GSM>/<GSM>.fastq.gz            (single-end)
#   <OUTDIR>/<GSM>/<GSM>_1.fastq.gz, _2.fastq.gz   (paired-end)
#
# Usage:
#   Sequential (small test, e.g. first 2 rows):
#     ./download.sh MANIFEST.csv OUTDIR --limit 2
#   Full run, single machine, sequential:
#     ./download.sh manifest.csv /path/to/outdir
#   SLURM array (optional -- only if you have a Slurm cluster; see download.slurm):
#     sbatch --array=1-N%5 download.slurm manifest.csv /path/to/outdir
#     (SLURM_ARRAY_TASK_ID selects one manifest row/GSM per task)

set -euo pipefail
IFS=$'\n\t'

MANIFEST="${1:?Usage: $0 MANIFEST.csv OUTDIR [--limit N | --array-mode]}"
OUTDIR="${2:?Usage: $0 MANIFEST.csv OUTDIR [--limit N | --array-mode]}"
MODE="${3:-}"
LIMIT="${4:-0}"

mkdir -p "$OUTDIR" logs

# manifest columns: gsm,kind,srx,gse,runs,fastq_bytes,fastq_ftp_by_run
# skip header, keep 1-indexed line numbers matching SLURM_ARRAY_TASK_ID

download_row() {
    local gsm="$1" ftp_by_run="$2"
    local sample_dir="$OUTDIR/$gsm"
    mkdir -p "$sample_dir"

    if [ -z "$ftp_by_run" ]; then
        echo "[$gsm] SKIP: no fastq_ftp_by_run in manifest (no ENA size/file record found)"
        return
    fi

    IFS='|' read -ra run_groups <<< "$ftp_by_run"
    local n_runs=${#run_groups[@]}

    # first run determines layout (single-end vs paired-end); all runs for a GSM
    # are expected to share the same layout
    IFS=';' read -ra first_urls <<< "${run_groups[0]}"
    local n_reads=${#first_urls[@]}

    for r in $(seq 0 $((n_reads - 1))); do
        local suffix=""
        if [ "$n_reads" -gt 1 ]; then
            suffix="_$((r+1))"
        fi
        local out="$sample_dir/${gsm}${suffix}.fastq.gz"
        if [ -s "$out" ]; then
            echo "[$gsm] already present: $out"
            continue
        fi

        local part_files=()
        for g in "${!run_groups[@]}"; do
            IFS=';' read -ra urls <<< "${run_groups[$g]}"
            local url="https://${urls[$r]}"
            local part="$sample_dir/${gsm}.run${g}${suffix}.fastq.gz.part"
            if [ ! -s "$part" ]; then
                echo "[$gsm] downloading run $((g+1))/$n_runs -> $part"
                curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors -o "$part.partial" "$url" && mv "$part.partial" "$part"
            fi
            part_files+=("$part")
        done

        if [ "$n_runs" -gt 1 ]; then
            echo "[$gsm] concatenating $n_runs runs -> $out"
        fi
        cat "${part_files[@]}" > "$out.partial"
        mv "$out.partial" "$out"
        rm -f "${part_files[@]}"
    done
}

if [ "$MODE" = "--array-mode" ]; then
    : "${SLURM_ARRAY_TASK_ID:?--array-mode requires SLURM_ARRAY_TASK_ID (run via sbatch --array=...)}"
    line=$((SLURM_ARRAY_TASK_ID + 1))  # +1 to skip header
    row=$(sed -n "${line}p" "$MANIFEST")
    if [ -z "$row" ]; then
        echo "No manifest row for array index $SLURM_ARRAY_TASK_ID (line $line) -- past end of file"
        exit 0
    fi
    IFS=',' read -r gsm kind srx gse runs fastq_bytes fastq_ftp_by_run <<< "$row"
    download_row "$gsm" "$fastq_ftp_by_run"
else
    n=0
    tail -n +2 "$MANIFEST" | while IFS=',' read -r gsm kind srx gse runs fastq_bytes fastq_ftp_by_run; do
        n=$((n+1))
        if [ "$MODE" = "--limit" ] && [ "$LIMIT" -gt 0 ] && [ "$n" -gt "$LIMIT" ]; then
            break
        fi
        download_row "$gsm" "$fastq_ftp_by_run"
    done
fi
