"""Assemble the final per-GSM download manifest from Step 7's run-level output.

One row per GSM, not per run: a GSM with more than one SRA run lists all of
them, grouped, rather than getting one manifest row per run. An earlier
version of this tool used one-row-per-run, which meant two downloader tasks
for the same multi-run GSM would race on the same output filename -- one
would fail outright, and worse, a "successful" task could silently overwrite
the other run's data, quietly losing half the reads. Grouping by GSM here and
having the downloader concatenate runs in order (see download.sh) fixes both.

Schema: gsm,kind,srx,gse,runs,fastq_bytes,fastq_ftp_by_run
  - runs:            ';'-joined SRA run accessions for this GSM
  - fastq_ftp_by_run: '|'-joined per-run URL groups; within a run, files are
                       ';'-joined (paired-end R1;R2)
  - fastq_bytes:      total bytes summed across every run/file for this GSM
"""
import argparse
import csv
import sys
from collections import defaultdict


def build(run_rows):
    by_gsm = defaultdict(list)
    for r in run_rows:
        by_gsm[r["gsm"]].append(r)

    out_rows = []
    for gsm, rs in by_gsm.items():
        rs = [r for r in rs if r["fastq_ftp"]]
        if not rs:
            continue
        kind = rs[0]["kind"]
        srx = rs[0]["srx"]
        gse = rs[0].get("gse", "")
        runs = ";".join(r["run_accession"] for r in rs)
        ftp_by_run = "|".join(r["fastq_ftp"] for r in rs)
        total_bytes = sum(int(b) for r in rs for b in r["fastq_bytes"].split(";") if r["fastq_bytes"])
        out_rows.append({
            "gsm": gsm, "kind": kind, "srx": srx, "gse": gse,
            "runs": runs, "fastq_bytes": total_bytes, "fastq_ftp_by_run": ftp_by_run,
        })
    out_rows.sort(key=lambda r: (r["gse"], r["gsm"]))
    return out_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="infile", required=True, help="run-level CSV from resolve_sra.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gse-map", default=None,
                     help="optional CSV (gsm,gse) to attach a series accession when --in lacks one "
                          "(resolve_sra.py's output doesn't carry gse; pass your Step 6 table here)")
    args = ap.parse_args()

    with open(args.infile, newline="") as f:
        run_rows = list(csv.DictReader(f))

    if args.gse_map:
        with open(args.gse_map, newline="") as f:
            gse_of = {r["gsm"]: r["gse"] for r in csv.DictReader(f)}
        for r in run_rows:
            r["gse"] = gse_of.get(r["gsm"], "")
    else:
        for r in run_rows:
            r.setdefault("gse", "")

    out_rows = build(run_rows)
    total_bytes = sum(r["fastq_bytes"] for r in out_rows)
    multi_run = sum(1 for r in out_rows if ";" in r["runs"])
    print(f"{len(out_rows)} GSMs in manifest ({multi_run} multi-run), {total_bytes / 1e12:.2f} TB total", file=sys.stderr)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gsm", "kind", "srx", "gse", "runs", "fastq_bytes", "fastq_ftp_by_run"],
                            lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
