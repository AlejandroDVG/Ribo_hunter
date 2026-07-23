"""Step 7: resolve each verified GSM to its SRA experiment (SRX) accession via
the same esummary record type used for discovery, then query the EBI ENA
Portal API directly for exact run-level byte sizes and FASTQ URLs -- no NCBI
SRA Toolkit needed, since ENA re-serves every SRA read as ready-to-use
gzipped FASTQ and reports its exact size up front.
"""
import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import common


def resolve_srx(gsm):
    uid = common.esearch_single_uid(gsm)
    time.sleep(common.NCBI_DELAY_S)
    if not uid:
        return gsm, ""
    recs = common.esummary([uid])
    rec = recs.get(uid, {})
    for rel in rec.get("extrelations", []):
        if rel.get("relationtype") == "SRA":
            return gsm, rel.get("targetobject", "")
    return gsm, ""


def resolve_all_srx(gsms, max_workers=6, progress=True):
    """Sequential, not threaded: this hits NCBI's eutils endpoints, which enforce
    a strict ~3 req/sec unauthenticated rate limit -- concurrent requests here
    (unlike verify_species.py's GEO acc.cgi calls, a separate, more tolerant host)
    trigger HTTP 429s. `max_workers` is accepted for CLI compatibility but ignored."""
    srx_map = {}
    for i, gsm in enumerate(gsms):
        _, srx = resolve_srx(gsm)
        srx_map[gsm] = srx
        if progress and (i + 1) % 50 == 0:
            print(f"...resolved SRX for {i + 1}/{len(gsms)}", file=sys.stderr, flush=True)
    return srx_map


def fetch_ena_runs(srx):
    return common.ena_filereport(srx)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="infile", required=True, help="CSV with a 'gsm' column and a 'kind' column (e.g. ribo-seq/rna-seq)")
    ap.add_argument("--gsm-col", default="gsm")
    ap.add_argument("--kind-col", default="kind")
    ap.add_argument("--out", required=True, help="run-level CSV: gsm,kind,srx,run_accession,fastq_bytes,fastq_ftp")
    ap.add_argument("--max-workers", type=int, default=6)
    args = ap.parse_args()

    with open(args.infile, newline="") as f:
        rows = list(csv.DictReader(f))
    gsm_kind = {r[args.gsm_col]: r.get(args.kind_col, "") for r in rows if r.get(args.gsm_col)}
    gsms = list(gsm_kind.keys())

    print(f"resolving SRX for {len(gsms)} GSMs...", file=sys.stderr)
    srx_map = resolve_all_srx(gsms, max_workers=args.max_workers)
    n_resolved = sum(1 for v in srx_map.values() if v)
    print(f"SRX resolved for {n_resolved}/{len(gsms)}", file=sys.stderr)

    out_rows = []
    resolved_srx = [(g, s) for g, s in srx_map.items() if s]
    done = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futures = {ex.submit(fetch_ena_runs, srx): (gsm, srx) for gsm, srx in resolved_srx}
        for fut in as_completed(futures):
            gsm, srx = futures[fut]
            try:
                runs = fut.result()
            except Exception as e:  # noqa: BLE001
                print(f"  ENA lookup failed for {gsm}/{srx}: {e}", file=sys.stderr)
                runs = []
            if not runs:
                out_rows.append({"gsm": gsm, "kind": gsm_kind[gsm], "srx": srx,
                                  "run_accession": "", "fastq_bytes": "", "fastq_ftp": ""})
            else:
                for rec in runs:
                    out_rows.append({
                        "gsm": gsm, "kind": gsm_kind[gsm], "srx": srx,
                        "run_accession": rec.get("run_accession", ""),
                        "fastq_bytes": rec.get("fastq_bytes", ""),
                        "fastq_ftp": rec.get("fastq_ftp", ""),
                    })
            done += 1
            if done % 50 == 0:
                print(f"...ENA resolved {done}/{len(resolved_srx)}", file=sys.stderr, flush=True)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gsm", "kind", "srx", "run_accession", "fastq_bytes", "fastq_ftp"])
        w.writeheader()
        w.writerows(out_rows)

    total_bytes = sum(
        int(b) for r in out_rows if r["fastq_bytes"] for b in r["fastq_bytes"].split(";")
    )
    print(f"wrote {args.out} -- {len(out_rows)} run rows, {total_bytes / 1e12:.2f} TB total", file=sys.stderr)


if __name__ == "__main__":
    main()
