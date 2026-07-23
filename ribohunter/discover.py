"""Step 1-2: discovery search + per-sample keyword confirmation.

A GEO series title/abstract describes the whole submission, which is often a
SuperSeries bundling several assay types. A series-level keyword hit is not
evidence that the series actually contains ribosome-profiling data -- the
match has to be confirmed on each individual sample's own title.
"""
import argparse
import csv
import re
import sys
import time

from . import common

DEFAULT_TERMS = ["ribosome profiling", "ribo-seq", "riboseq"]

# Broader than the discovery terms on purpose: catches synonyms depositors use
# in sample titles even when the series-level abstract used different wording.
SAMPLE_TITLE_RE = re.compile(
    r"ribo[-\s]?seq|ribosome\s?profil\w*|ribosome[-\s]?protected\w*|footprint\w*"
    r"|(?<![A-Za-z])RPFs?(?![A-Za-z])|translatom\w*",
    re.I,
)


def build_discovery_term(organism, extra_exclude_terms=None):
    ors = " OR ".join(f'"{t}"[All Fields]' if " " in t or "-" in t else f"{t}[All Fields]" for t in DEFAULT_TERMS)
    term = f"({ors}) AND \"{organism}\"[Organism] AND \"gse\"[Filter]"
    return term


def discover_series(organism):
    """Step 1: find candidate GEO series (GSE UIDs) for an organism."""
    term = build_discovery_term(organism)
    uids = common.esearch(term)
    return uids


def confirm_samples(series_uids, organism, progress=True):
    """Step 2: fetch full series records and keep only samples whose OWN title
    matches ribosome-profiling terminology. Returns a list of dicts, one per
    confirmed sample, with the series accession, sample accession, and title."""
    confirmed = []
    n = len(series_uids)
    records = common.esummary(series_uids)
    for i, uid in enumerate(series_uids):
        rec = records.get(uid)
        if not rec:
            continue
        gse_acc = rec.get("accession", "")
        for sample in rec.get("samples", []):
            title = sample.get("title", "")
            gsm_acc = sample.get("accession", "")
            if not gsm_acc or not title:
                continue
            if SAMPLE_TITLE_RE.search(title):
                confirmed.append({"gse": gse_acc, "gsm": gsm_acc, "title": title})
        if progress and (i + 1) % 100 == 0:
            print(f"...confirmed {i + 1}/{n} series", file=sys.stderr, flush=True)
    return confirmed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, help='GEO organism string, e.g. "Drosophila melanogaster"')
    ap.add_argument("--out", required=True, help="output CSV: gse,gsm,title")
    args = ap.parse_args()

    print(f"Step 1: discovery search for organism={args.organism!r}", file=sys.stderr)
    series_uids = discover_series(args.organism)
    print(f"  {len(series_uids)} candidate series", file=sys.stderr)

    print("Step 2: per-sample title confirmation", file=sys.stderr)
    confirmed = confirm_samples(series_uids, args.organism)
    confirmed_series = {c["gse"] for c in confirmed}
    print(f"  {len(confirmed)} confirmed samples across {len(confirmed_series)} series", file=sys.stderr)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gse", "gsm", "title"])
        w.writeheader()
        w.writerows(confirmed)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
