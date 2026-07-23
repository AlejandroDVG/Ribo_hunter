"""Step 3: deduplicate cross-listed samples.

GEO SuperSeries and SubSeries can both list the same physical sample
accession. Group by GSM; where one appears under more than one GSE, keep a
single row and record the alternate accession(s) instead of counting the
sample twice.
"""
import argparse
import csv
import sys
from collections import OrderedDict


def dedup(rows):
    """rows: list of dicts with at least 'gsm' and 'gse'. Returns deduplicated
    rows (one per unique gsm, first-seen gse kept as primary) each with an
    added 'also_listed_under' field (semicolon-joined alternate GSEs, or '')."""
    by_gsm = OrderedDict()
    for row in rows:
        gsm = row["gsm"]
        by_gsm.setdefault(gsm, []).append(row)

    out = []
    for gsm, entries in by_gsm.items():
        primary = entries[0]
        alt_gses = [e["gse"] for e in entries[1:] if e["gse"] != primary["gse"]]
        primary = dict(primary)
        primary["also_listed_under"] = ";".join(dict.fromkeys(alt_gses))
        out.append(primary)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="infile", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.infile, newline="") as f:
        rows = list(csv.DictReader(f))

    deduped = dedup(rows)
    n_dupe_gsms = sum(1 for r in deduped if r["also_listed_under"])
    print(f"{len(rows)} rows -> {len(deduped)} unique GSMs ({n_dupe_gsms} cross-listed under >1 series)", file=sys.stderr)

    fieldnames = list(rows[0].keys()) + ["also_listed_under"] if rows else ["gse", "gsm", "title", "also_listed_under"]
    # preserve original column order, add the new one at the end without duplicating it
    fieldnames = list(dict.fromkeys(fieldnames))
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(deduped)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
