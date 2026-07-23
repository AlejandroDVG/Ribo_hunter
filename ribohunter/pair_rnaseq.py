"""Step 6: pair each ribosome-profiling sample with its matched RNA-seq control.

Match by normalizing both titles -- stripping assay-type tokens -- and
requiring an exact match on what's left, within the same series (a sample's
own GSE plus anything it's cross-listed under from Step 3's dedup).

This step needs its OWN discovery pass over each series' full sample list
(not just the already-confirmed ribo-seq samples), because the RNA-seq
partner's title generally does NOT match the ribo-seq keyword regex at all --
that's the whole point of it being the RNA-seq control.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

from . import common

# Tokens stripped before matching. Asymmetric stripping (e.g. recognizing
# "footprint" but not "protected RNA", or missing bare "Ribo" outside
# "Ribo-seq") is exactly the bug class that caused missed pairs in prior runs
# -- keep this list as wide as plausible rather than adding tokens reactively.
STRIP_TOKENS = [
    # "seq" is almost always joined to "ribo"/"rna" with a hyphen/space/underscore
    # ("Ribo-seq", "RNA-seq") rather than directly ("riboseq") -- \w* doesn't cross
    # that separator, so an unbridged pattern here silently fails to strip the far
    # more common hyphenated form and breaks pairing for nearly every real title.
    r"ribo[\s_-]*seq", r"m?rna[\s_-]*seq", r"\binput\b", r"protected\s*rna", r"footprint\w*",
    r"\bribo\b", r"\bm?rna\b", r"rpfs?\b", r"mono[-\s]?somes?", r"di[-\s]?somes?", r"total\s*rna",
    r"\brep(licate)?\s*\d*\b", r"\brun\s*\d*\b",
]
STRIP_RE = re.compile("|".join(STRIP_TOKENS), re.I)
REPLICATE_RE = re.compile(r"rep(?:licate)?\s*(\d+)", re.I)


def normalize_title(title):
    t = re.sub(r"\[r[_-]", "[", title, flags=re.I)  # e.g. "[r_pR1]" vs sibling "[pR1]" bracket-naming mismatch
    t = STRIP_RE.sub(" ", t)
    t = re.sub(r"[^a-z0-9]+", " ", t.lower())
    return t.strip()


def replicate_of(title):
    m = REPLICATE_RE.search(title)
    return m.group(1) if m else None


def fetch_series_samples(gse):
    """All samples (accession + title) belonging to a GSE, via esummary."""
    uids = common.esearch(f"{gse}[ACCN] AND GSE[Filter]")
    if not uids:
        return []
    recs = common.esummary(uids)
    for rec in recs.values():
        if rec.get("accession") == gse:
            return rec.get("samples", [])
    return []


def pair_series(ribo_rows, all_series_samples):
    """ribo_rows: confirmed ribo-seq rows (with 'gsm','title') for one series.
    all_series_samples: full sample list for that series (accession/title dicts).
    Returns dict: ribo gsm -> (paired gsm, paired title) or (None, None)."""
    ribo_gsms = {r["gsm"] for r in ribo_rows}
    candidates = [s for s in all_series_samples if s["accession"] not in ribo_gsms]
    by_norm = defaultdict(list)
    for s in candidates:
        by_norm[normalize_title(s["title"])].append(s)

    claimed = set()
    pairs = {}
    for r in ribo_rows:
        norm = normalize_title(r["title"])
        options = [s for s in by_norm.get(norm, []) if s["accession"] not in claimed]
        if not options:
            pairs[r["gsm"]] = (None, None)
            continue
        # precision guard: prefer a replicate-number match when there's more than one option left
        rep = replicate_of(r["title"])
        if rep and len(options) > 1:
            options = sorted(options, key=lambda s: replicate_of(s["title"]) != rep)
        chosen = options[0]
        claimed.add(chosen["accession"])
        pairs[r["gsm"]] = (chosen["accession"], chosen["title"])
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="infile", required=True, help="CSV with gsm,gse,title (+ also_listed_under)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.infile, newline="") as f:
        rows = list(csv.DictReader(f))

    by_series = defaultdict(list)
    for r in rows:
        by_series[r["gse"]].append(r)

    all_pairs = {}
    n_series = len(by_series)
    for i, (gse, ribo_rows) in enumerate(by_series.items()):
        samples = fetch_series_samples(gse)
        pairs = pair_series(ribo_rows, samples)
        all_pairs.update(pairs)
        if (i + 1) % 20 == 0:
            print(f"...paired {i + 1}/{n_series} series", file=sys.stderr, flush=True)

    n_paired = sum(1 for v in all_pairs.values() if v[0])
    print(f"paired {n_paired}/{len(rows)} ribo-seq samples ({n_paired / len(rows):.0%})", file=sys.stderr)

    fieldnames = list(rows[0].keys()) + ["Paired_RNAseq_sample_id", "Paired_RNAseq_sample_name"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            gsm, title = all_pairs.get(r["gsm"], (None, None))
            out_row = dict(r)
            out_row["Paired_RNAseq_sample_id"] = gsm or ""
            out_row["Paired_RNAseq_sample_name"] = title or ""
            w.writerow(out_row)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
