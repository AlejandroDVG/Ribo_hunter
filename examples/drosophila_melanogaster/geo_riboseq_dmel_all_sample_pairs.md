# Drosophila melanogaster Ribo-seq discovery — worked example

Run through `ribohunter`'s full 8-step method (see the top-level
[README](../../README.md)) with `--organism "Drosophila melanogaster"`.

## Counts at each step

| Step | Count |
|---|---:|
| 1. Candidate series (discovery search) | 22 |
| 2. Series with ≥1 sample confirmed by title | 10 |
| 2. Confirmed samples (pre-dedup) | 95 |
| 3. Unique GSMs after deduplication | 95 (0 cross-listed under >1 series) |
| 4. Confirmed *Drosophila melanogaster* | 79 |
| 4. Excluded — wrong organism | 16 |
| 6. Paired with an RNA-seq control | 30 / 79 (38%) |
| 5. Treatment annotated | 39 / 79 (49%) |
| 5. Adapter sequence / kit only | 0 / 40 |
| 7-8. Manifest size (ribo + paired rna, ENA-resolved) | 100 GSMs, 0.19 TB |

## Bug caught by step 4 (species re-verification)

Of the 95 samples that matched ribosome-profiling keywords at the per-sample
title level, **16 turned out not to be *Drosophila melanogaster*** when their
own GEO record was checked directly (see
[`species_contamination_excluded.csv`](species_contamination_excluded.csv)):

| Sample | Actual organism |
|---|---|
| GSM1276538 | Schizosaccharomyces pombe |
| GSM1276540, GSM1276544 | Mus musculus |
| GSM1276542 | Homo sapiens |
| GSM1276554–GSM1276562 (9 samples) | Danio rerio |
| GSM1276564, GSM1276566, GSM1316826 | Xenopus laevis |

All 16 came from a single multi-species comparative series tagged
`Drosophila melanogaster` at the series level — the series itself compares
translation across five organisms, and the per-series organism tag doesn't
reflect that. This is exactly the failure mode the README describes: had this
tool trusted the series-level tag (as a naive single-query search would),
all 16 would have been silently included as fly samples.

Two more samples matched the ribosome-profiling keyword regex only because
their title names the parent *study* ("... input mRNA for Ribosome
Profiling") while the sample itself is the RNA-seq input control, not a
ribo-seq assay — these were excluded from the confirmed set for the same
reason: a real ribo-specific term (ribo-seq, RPF, footprint,
ribosome-protected) has to appear, not just a mention of the study's name.

## Pairing

30 of 79 confirmed Ribo-seq samples (38%) were paired with a matched RNA-seq
control. The remaining 49 were checked individually rather than assumed
unpairable:

- Several series (e.g. `GSE166408`) have no RNA-seq samples in the series at
  all — genuinely unpaired, not a matching failure.
- Some series pair Ribo-seq "IP" samples against "Input" samples from the
  *same* ribosome-profiling protocol rather than against a separate RNA-seq
  assay (e.g. `GSE245380`'s "Ribo-Seq Ctl IP"/"Ribo-Seq Ctl Input" naming) —
  there's no separate RNA-seq sample to pair against.
- A couple of series number Ribo-seq and RNA-seq replicates with entirely
  different naming schemes and no shared identifier (e.g. `GSE153346`'s
  numeric sample IDs, `GSE52799`'s descriptive fraction names) — left
  unpaired rather than guessed at.

One systematic bug was found and fixed while building this pairing: the
token-stripping regex used `ribo\w*seq`/`rna\w*seq`, which doesn't cross the
hyphen in "Ribo-seq"/"RNA-seq" — the far more common way titles actually
write it (`mRNA-seq of S2 cells...` vs `Ribo-seq of S2 cells...` normalized
to `mrna seq of...` vs `seq of...` — not equal). Fixed to
`ribo[\s_-]*seq`/`m?rna[\s_-]*seq`, plus a bare `\bm?rna\b` strip for titles
using "mRNA" without "-seq" attached (e.g. `mRNA 0-2 hour embryo replicate A`
vs `footprint 0-2 hour embryo replicate A`). This raised pairing from 28 to
30 matched pairs and is now the default behavior, not something specific to
this example.

## Files in this directory

| File | Description |
|---|---|
| `geo_riboseq_dmel_all_sample_pairs.csv` | Final per-sample table — 79 confirmed Ribo-seq samples, steps 1-6 output. |
| `download_manifest_dmel.csv` | One row per GSM (ribo + paired rna-seq), ENA byte sizes and FASTQ URLs — step 7-8 output. |
| `species_contamination_excluded.csv` | The 16 samples dropped by step 4, with their actual (non-fly) organism. |

## Reproducing this

```bash
python3 -m ribohunter.discover --organism "Drosophila melanogaster" --out confirmed.csv
python3 -m ribohunter.dedup --in confirmed.csv --out deduped.csv
python3 -m ribohunter.verify_species --in deduped.csv --organism "Drosophila melanogaster" \
    --out verified.csv --excluded-out excluded.csv
python3 -m ribohunter.extract_metadata --in verified.csv --out metadata.csv
python3 -m ribohunter.pair_rnaseq --in metadata.csv --out paired.csv
# then build a gsm,kind,gse input covering both ribo samples and their paired rna-seq
# samples, resolve_sra.py + build_manifest.py as shown in the top-level README
```
