"""Step 4: per-sample species re-verification.

The organism tag used for discovery is set at the SERIES level and is not
reliable at the sample level -- SuperSeries frequently bundle more than one
organism. Fetch each sample's own full-text GEO record and read
Sample_organism_ch1 directly rather than trusting the series-level tag.
This step has caught real cross-contamination in every organism run so far
(both directions: target-organism samples missing it, and other-organism
samples slipping in) -- do not skip it.
"""
import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import common


def fetch_organism(gsm):
    try:
        fields = common.fetch_geo_full_record(gsm)
    except Exception as e:  # noqa: BLE001
        return gsm, {"organism": "", "error": str(e)}
    organism = fields.get("Sample_organism_ch1", [""])[0]
    return gsm, {"organism": organism, "error": ""}


def verify_all(gsms, organism, max_workers=6, progress=True):
    """Returns (confirmed, excluded) -- each a list of dicts with gsm/organism."""
    confirmed, excluded = [], []
    n = len(gsms)
    done = 0
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(fetch_organism, gsm): gsm for gsm in gsms}
        for fut in as_completed(futures):
            gsm, info = fut.result()
            row = {"gsm": gsm, "organism": info["organism"]}
            if info["organism"].strip().lower() == organism.strip().lower():
                confirmed.append(row)
            else:
                excluded.append(row)
            done += 1
            if progress and done % 50 == 0:
                print(f"...verified {done}/{n}", file=sys.stderr, flush=True)
    return confirmed, excluded


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="infile", required=True, help="CSV with a 'gsm' column")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True, help="output CSV: only samples confirmed as --organism")
    ap.add_argument("--excluded-out", default=None, help="optional CSV of samples dropped for wrong organism")
    ap.add_argument("--max-workers", type=int, default=6)
    args = ap.parse_args()

    with open(args.infile, newline="") as f:
        rows = list(csv.DictReader(f))

    gsms = [r["gsm"] for r in rows]
    print(f"verifying organism for {len(gsms)} samples against GEO's own per-sample record...", file=sys.stderr)
    confirmed, excluded = verify_all(gsms, args.organism, max_workers=args.max_workers)
    print(f"confirmed {args.organism}: {len(confirmed)} | excluded (other organism): {len(excluded)}", file=sys.stderr)
    if excluded:
        print("excluded samples and their actual organism:", file=sys.stderr)
        for row in excluded:
            print(f"  {row['gsm']}: {row['organism'] or 'FETCH_ERROR'}", file=sys.stderr)

    confirmed_gsms = {r["gsm"] for r in confirmed}
    out_rows = [r for r in rows if r["gsm"] in confirmed_gsms]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {args.out}", file=sys.stderr)

    if args.excluded_out:
        with open(args.excluded_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["gsm", "organism"])
            w.writeheader()
            w.writerows(excluded)
        print(f"wrote {args.excluded_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
