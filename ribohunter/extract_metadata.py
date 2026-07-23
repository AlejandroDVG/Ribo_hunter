"""Step 5: metadata extraction (cell line/tissue, treatment, adapter) from each
sample's own full-text GEO record (Sample_characteristics_ch1,
Sample_extract_protocol_ch1, Sample_data_processing, Sample_title).

Automatic extraction here is deliberately conservative: it surfaces GEO's own
structured fields and canonicalizes obvious matches (known drug names, an
explicit adapter sequence), but it does NOT attempt full manual homogenization
of cell-line/tissue naming across a whole corpus (e.g. collapsing "HeLa cells",
"Hela Cells", "HELA" variants into one canonical label) -- that alias table is
corpus- and organism-specific and was built by hand for the human/mouse runs
this tool was extracted from. Treat Cell_line_source_raw as a starting point
for your own cleanup pass, not a finished column.
"""
import argparse
import csv
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import common

# ---------- treatment ----------

DRUG_PATTERNS = [
    (re.compile(r"harringtonine|\bHHT\b", re.I), "Harringtonine"),
    (re.compile(r"lactimidomycin|\bLTM\b", re.I), "Lactimidomycin (LTM)"),
    (re.compile(r"cycloheximide|cyclohexamide|\bCHX\b", re.I), "Cycloheximide (CHX)"),
    (re.compile(r"puromycin", re.I), "Puromycin"),
    (re.compile(r"anisomycin", re.I), "Anisomycin"),
    (re.compile(r"emetine", re.I), "Emetine"),
    (re.compile(r"actinomycin\s*D", re.I), "Actinomycin D"),
    (re.compile(r"thapsigargin", re.I), "Thapsigargin"),
    (re.compile(r"tunicamycin", re.I), "Tunicamycin"),
    (re.compile(r"sodium arsenite|\bNaAs\b|arsenite", re.I), "Sodium arsenite (NaAsO2)"),
    (re.compile(r"\bdoxycycline\b|\bdox\b(?!orubicin)", re.I), "Doxycycline"),
    (re.compile(r"\bDMSO\b", re.I), "DMSO (vehicle control)"),
    (re.compile(r"\brapamycin\b", re.I), "Rapamycin"),
    (re.compile(r"\btorin", re.I), "Torin"),
    (re.compile(r"\bLPS\b", re.I), "LPS"),
    (re.compile(r"\btamoxifen\b", re.I), "Tamoxifen"),
    (re.compile(r"interferon|\bIFN", re.I), "Interferon"),
    (re.compile(r"\bcontrol\b|\buntreated\b|\bmock\b|\bvehicle\b", re.I), "Control/untreated"),
]
CONFIRMED_NONE = {"none", "na", "n/a", "n.a.", "not applicable", "no drug", "nan"}
NO_SIGNAL = {"unknown", ""}
TREATMENT_KEY_PATTERN = re.compile(r"treat|^agent$|sirna", re.I)
VALUE_TREATMENT_HINT = re.compile(r"\btreatment\b|\btreated\b", re.I)
# withheld, not administered -- "- doxycycline" / "without dox" mean the opposite of a match on the drug name
NEGATION_PATTERN = re.compile(r"^\s*-\s*\w|\(\s*-\s*\)|\bwithout\b", re.I)


def _parse_treatment_values(characteristics):
    vals = []
    for part in characteristics.split("|"):
        if ":" not in part:
            continue
        key, _, val = part.partition(":")
        key, val = key.strip().lower(), val.strip()
        if not val:
            continue
        if TREATMENT_KEY_PATTERN.search(key):
            vals.append(val)
        elif any(pat.search(key) for pat, label in DRUG_PATTERNS if label != "Control/untreated"):
            # the drug name itself used as the characteristic key, e.g. "cycloheximide: 0.1 mg/mL"
            vals.append(f"{key} {val}")
        elif VALUE_TREATMENT_HINT.search(val):
            vals.append(val)
    return vals


def _canonicalize_treatment(val):
    low = val.strip().lower()
    if low in CONFIRMED_NONE:
        return ["Control/untreated"]
    if low in NO_SIGNAL:
        return []
    if NEGATION_PATTERN.search(val):
        return ["Control/untreated"]
    matches = list(dict.fromkeys(label for pat, label in DRUG_PATTERNS if pat.search(val)))
    return matches or [val]


def extract_treatment(characteristics, title):
    raw_values = _parse_treatment_values(characteristics)
    out = []
    for v in raw_values:
        for c in _canonicalize_treatment(v):
            if c not in out:
                out.append(c)
    if out:
        return "; ".join(out)
    for pat, label in DRUG_PATTERNS:
        if pat.search(title):
            return label
    return ""


# ---------- adapter ----------

SEQ_NEAR_PAREN = re.compile(r"adapt[oe]r[^()]{0,60}\(([ACGTUacgtu]{8,60})\)", re.I)
SEQ_AFTER_COLON = re.compile(r"adapt[oe]r[^:]{0,40}:\s*([ACGTUacgtu]{8,60})\b", re.I)
SEQ_LOOSE = re.compile(r"adapt[oe]r\s+sequence[^A-Za-z]{0,20}([ACGTUacgtu]{8,60})\b", re.I)
KIT_PATTERNS = [
    re.compile(r"TruSeq\s+Ribo\s*Profile", re.I),
    re.compile(r"TruSeq\s+Small\s*RNA", re.I),
    re.compile(r"NEBNext\s+Small\s*RNA", re.I),
    re.compile(r"NEBNext\s+Multiplex\s+Small\s*RNA", re.I),
    re.compile(r"Illumina\s+TruSeq", re.I),
    re.compile(r"CleanTag", re.I),
    re.compile(r"NextFlex\s+Small\s*RNA", re.I),
]


def extract_adapter(text):
    seq = ""
    for pat in (SEQ_NEAR_PAREN, SEQ_AFTER_COLON, SEQ_LOOSE):
        m = pat.search(text)
        if m:
            candidate = m.group(1).upper().replace("U", "T")
            if len(set(candidate)) > 1:  # not just a homopolymer run misfiring the regex
                seq = candidate
                break
    kit = ""
    for pat in KIT_PATTERNS:
        m = pat.search(text)
        if m:
            kit = m.group(0)
            break
    return seq, kit


# ---------- cell line / tissue (raw only, see module docstring) ----------

CELL_LINE_KEY_RE = re.compile(r"cell[\s_-]?line|tissue|cell[\s_-]?type", re.I)


def extract_cell_line_raw(characteristics, source_name):
    for part in characteristics.split("|"):
        if ":" not in part:
            continue
        key, _, val = part.partition(":")
        if CELL_LINE_KEY_RE.search(key.strip()) and val.strip():
            return val.strip()
    return source_name.strip()


# ---------- driver ----------


def _clean(s):
    return re.sub(r"\s+", " ", s.replace("\r", "").replace("\xa0", " ")).strip()


def fetch_one(gsm):
    fields = common.fetch_geo_full_record(gsm)
    title = _clean(fields.get("Sample_title", [""])[0])
    characteristics = _clean(" | ".join(fields.get("Sample_characteristics_ch1", [])))
    extract_protocol = _clean(" | ".join(fields.get("Sample_extract_protocol_ch1", [])))
    data_processing = _clean(" | ".join(fields.get("Sample_data_processing", [])))
    source_name = _clean(fields.get("Sample_source_name_ch1", [""])[0])

    combined_protocol_text = extract_protocol + " | " + data_processing
    seq, kit = extract_adapter(combined_protocol_text)
    return {
        "gsm": gsm,
        "title": title,
        "Cell_line_source": extract_cell_line_raw(characteristics, source_name),
        "Cell_line_source_raw": source_name,
        "Treatment": extract_treatment(characteristics, title),
        "Adapter_sequence": seq,
        "Adapter_kit_mentioned": kit,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="infile", required=True, help="CSV with a 'gsm' column")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-workers", type=int, default=6)
    args = ap.parse_args()

    with open(args.infile, newline="") as f:
        rows = list(csv.DictReader(f))
    gsms = [r["gsm"] for r in rows]

    print(f"extracting metadata for {len(gsms)} samples...", file=sys.stderr)
    results = {}
    done = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futures = {ex.submit(fetch_one, gsm): gsm for gsm in gsms}
        for fut in as_completed(futures):
            gsm = futures[fut]
            try:
                results[gsm] = fut.result()
            except Exception as e:  # noqa: BLE001
                print(f"  fetch failed for {gsm}: {e}", file=sys.stderr)
                results[gsm] = {"gsm": gsm, "title": "", "Cell_line_source": "", "Cell_line_source_raw": "",
                                 "Treatment": "", "Adapter_sequence": "", "Adapter_kit_mentioned": ""}
            done += 1
            if done % 50 == 0:
                print(f"...{done}/{len(gsms)}", file=sys.stderr, flush=True)

    n_treat = sum(1 for r in results.values() if r["Treatment"])
    n_seq = sum(1 for r in results.values() if r["Adapter_sequence"])
    n_kit = sum(1 for r in results.values() if r["Adapter_kit_mentioned"])
    print(f"treatment annotated: {n_treat}/{len(gsms)} | adapter seq: {n_seq} | kit only: {n_kit}", file=sys.stderr)

    meta_fields = ["Cell_line_source", "Cell_line_source_raw", "Treatment", "Adapter_sequence", "Adapter_kit_mentioned"]
    fieldnames = list(rows[0].keys()) + meta_fields
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            meta = results.get(row["gsm"], {})
            out_row = dict(row)
            for mf in meta_fields:
                out_row[mf] = meta.get(mf, "")
            w.writerow(out_row)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
