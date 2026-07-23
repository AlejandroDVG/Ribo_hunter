# ribohunter

Find, verify, and download every public ribosome-profiling (Ribo-seq) dataset
on GEO for a given organism.

`ribohunter` automates the eight-step method below end to end: search GEO,
confirm each hit at the *sample* level (not just the series), re-verify every
sample's organism against its own record, extract what metadata is available
(cell line/tissue, treatment, adapter), pair each Ribo-seq sample with its
matched RNA-seq control where one exists, resolve everything to exact
sequencing-run sizes and URLs via ENA, and build a manifest + downloader.

The worked example in this repo is *Drosophila melanogaster* — see
[`examples/drosophila_melanogaster/`](examples/drosophila_melanogaster/).

## Why not just search GEO once and download the hits?

Because a single series-level search silently gets both false positives and
false negatives:

- **GEO SuperSeries bundle multiple assay types and organisms.** A series
  whose abstract mentions "ribosome profiling" may be mostly RNA-seq with one
  Ribo-seq subseries — or worse, a multi-species comparative study where only
  some samples are your target organism. Running this tool's species
  re-verification step on the Drosophila example found **16 of 95
  keyword-matched samples were not actually *Drosophila melanogaster*** —
  Xenopus laevis, Danio rerio, Mus musculus, Schizosaccharomyces pombe, and
  Homo sapiens samples, all pulled in via a multi-organism comparative series
  that matched at the series level. This is not a one-off: the same check
  found analogous contamination (in both directions) on every organism this
  tool has been run against so far.
- **The organism tag used for discovery is set at the series level and just
  isn't reliable per-sample.** You have to fetch each sample's own record.
- **Series titles/abstracts don't always name the assay.** Some series have
  no ribosome-profiling language anywhere except in individual sample titles.

`ribohunter` exists because getting this right requires the same handful of
gotchas to be handled every time, for every organism — so it's a tool
instead of a one-off script per corpus.

## Install

```
git clone <this repo>
cd ribohunter
pip install -r requirements.txt   # just pandas; everything else is stdlib
```

No API key is required for the request volume this tool uses (NCBI allows
~3 unauthenticated requests/second; the tool throttles itself to stay under
that). No SRA Toolkit is required either — ENA re-serves every SRA-deposited
read as a ready-to-use gzipped FASTQ file and reports its exact size up
front, so downloads go straight from ENA's FTP.

## Usage

Each step is its own module and writes a CSV that the next step reads. Run
them in order for a new organism:

```bash
ORGANISM="Drosophila melanogaster"
OUT=examples/drosophila_melanogaster   # pick your own output dir

# 1-2: discovery search + per-sample keyword confirmation
python3 -m ribohunter.discover --organism "$ORGANISM" --out $OUT/confirmed.csv

# 3: deduplicate cross-listed samples (SuperSeries/SubSeries)
python3 -m ribohunter.dedup --in $OUT/confirmed.csv --out $OUT/deduped.csv

# 4: per-sample species re-verification (do not skip this)
python3 -m ribohunter.verify_species --in $OUT/deduped.csv --organism "$ORGANISM" \
    --out $OUT/verified.csv --excluded-out $OUT/excluded.csv

# 5: metadata extraction (cell line/tissue, treatment, adapter)
python3 -m ribohunter.extract_metadata --in $OUT/verified.csv --out $OUT/metadata.csv

# 6: pair each ribo-seq sample with its RNA-seq control, where one exists
python3 -m ribohunter.pair_rnaseq --in $OUT/metadata.csv --out $OUT/paired.csv

# 7: resolve every sample (ribo + paired rna) to SRA/ENA run-level sizes and URLs
#    (build your own gsm,kind,gse input -- one row per ribo sample plus one row
#    per paired rna-seq sample; see examples/ for the exact shape)
python3 -m ribohunter.resolve_sra --in $OUT/all_gsms.csv --out $OUT/runs.csv

# assemble the final per-GSM download manifest (multi-run-safe, see below)
python3 -m ribohunter.build_manifest --in $OUT/runs.csv --gse-map $OUT/gse_map.csv \
    --out $OUT/download_manifest.csv

# 8: download -- sequential by default, or as an optional Slurm array
mkdir -p logs
bash ribohunter/download.sh $OUT/download_manifest.csv /path/to/fastq_out --limit 2   # test first
bash ribohunter/download.sh $OUT/download_manifest.csv /path/to/fastq_out             # full run
# optional, only if you have a Slurm cluster:
sbatch --array=1-N%5 ribohunter/download.slurm $OUT/download_manifest.csv /path/to/fastq_out
```

## The 8 steps, and what actually goes wrong at each one

### 1-2. Discovery + per-sample keyword confirmation

Query GEO (`esearch`, database `gds`) for series matching ribosome-profiling
terminology, restricted to the target organism and to series-type records:

```
(ribosome profiling[All Fields] OR "ribo-seq"[All Fields] OR riboseq[All Fields])
AND "<organism>"[Organism] AND "gse"[Filter]
```

Keep the term list short and generic — recall matters more here than
precision, because the next step does the precision work. A series-level
keyword hit is not evidence the series contains Ribo-seq data: SuperSeries
routinely bundle an RNA-seq-only subseries under an abstract that only
mentions a *sibling* Ribo-seq subseries. `ribohunter` re-fetches full series
records and requires the keyword match on **each individual sample title**:

```
ribo[-\s]?seq | ribosome\s?profil\w* | ribosome[-\s]?protected\w* | footprint\w* | RPF | translatom\w*
```

Some series have no ribosome-profiling language anywhere in the series title
at all — the assay is only named in individual sample titles. A series-title-
only filter would silently miss these.

### 3. Deduplicate cross-listed samples

GEO SuperSeries and SubSeries can both list the same physical sample
accession. Group by GSM; where one appears under more than one GSE, keep a
single row and record the alternate accession(s) instead of counting the
sample twice.

### 4. Species re-verification, per sample — do not skip this

The organism tag used for discovery is set at the *series* level, not the
sample level. Fetch each sample's own full-text GEO record and read
`Sample_organism_ch1` directly:

```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSM}&targ=self&form=text&view=full
```

This has caught real contamination on every organism run so far, in both
directions: target-organism samples that the discovery search missed because
the series was tagged with a different organism, and other-organism samples
that slipped in because the series was tagged as the target organism. See the
Drosophila example above for a concrete case.

### 5. Metadata extraction: cell line/tissue, treatment, adapter

Extracted from the same per-sample record fetched in step 4
(`Sample_characteristics_ch1`, `Sample_extract_protocol_ch1`,
`Sample_data_processing`) because depositors write this information
inconsistently and no single field covers every case.

- **Cell line/tissue**: GEO's own structured `cell line:`/`tissue:`
  characteristic is preferred over the free-text source-name field when both
  exist. `ribohunter` extracts this raw, per-sample value but does **not**
  attempt full manual homogenization across a corpus (collapsing "HeLa
  cells"/"Hela Cells"/"HELA" into one canonical label) — that's a corpus- and
  organism-specific cleanup pass worth doing by hand once you have the raw
  values, not something this tool guesses at automatically.
- **Treatment**: matches characteristic keys containing `treat`/`agent`/
  `sirna`, canonicalizes against a known drug list, and keeps every drug
  present in a combined value rather than just the first match. Values like
  `"- doxycycline"` or `"without dox"` mean the drug was *withheld* — a naive
  keyword match would report the opposite of what happened, so these are
  explicitly caught and relabeled `Control/untreated`.
- **Adapter sequence**: looks for an explicit nucleotide sequence near the
  word "adapter" in the protocol/processing text. Where none is given, a
  recognized kit name is recorded as a separate field — a kit name is never
  converted into a guessed sequence, since kit-to-adapter mapping isn't fixed
  across kit versions and a wrong guess here silently corrupts every
  downstream trimming step.

### 6. Pair each sample with its matched RNA-seq control

Matches a Ribo-seq sample to its RNA-seq counterpart by normalizing both
titles — stripping assay-type tokens (`Ribo-seq`, `RNA-seq`/`mRNA-seq`,
`input`, `protected RNA`, bare `Ribo`/`RNA`, `footprint`, `monosome`/
`disome`, replicate/run numbers) — and requiring an exact match on what's
left, within the same series.

The tokens have to be stripped in a way that tolerates however the title
actually punctuates them: an early version of this stripping only matched
`ribo\w*seq`/`rna\w*seq`, which doesn't cross the hyphen in "Ribo-seq"/
"RNA-seq" — the far more common way people actually write it — and silently
failed to pair almost everything until fixed. If you're adapting this
pairing logic, audit it both ways before trusting it: check that no RNA-seq
sample gets claimed by more than one Ribo-seq sample (precision), and check
every unpaired Ribo-seq sample's series for an RNA-seq candidate that was
never claimed (recall). Not every unpaired sample is a bug — some series
genuinely have no RNA-seq counterpart, or number replicates with no shared
identifier at all.

### 7. Resolve to SRA/ENA and compute exact size

Each verified GSM resolves to an SRA experiment accession via `esummary`'s
`extrelations` (relation type `SRA`). Query the ENA Portal API directly with
that accession — no SRA Toolkit needed:

```
GET https://www.ebi.ac.uk/ena/portal/api/filereport
    ?accession=SRX...&result=read_run
    &fields=run_accession,fastq_bytes,fastq_ftp&format=json
```

This step hits NCBI's `esearch`/`esummary` endpoints, which enforce a strict
unauthenticated rate limit (~3 req/sec) — `resolve_sra.py` resolves SRX
accessions sequentially with an explicit delay rather than concurrently, to
avoid HTTP 429s (the earlier per-sample GEO-record fetches in step 4/5 hit a
different, more tolerant host and can safely run with a small thread pool).

### 8. Download, multi-run-safe

`build_manifest.py` assembles one manifest row **per GSM**, not per
sequencing run. A GSM with more than one SRA run lists all of them together;
`download.sh` downloads each run to its own temp file and concatenates them
in order into the final file. This matters: an earlier version of this
pipeline used one manifest row per run, so two downloader tasks for the same
multi-run GSM would race on the same output filename — one would fail
outright, and a "successful" task could silently overwrite the other run's
data, quietly losing half the reads with no error at all. Grouping by GSM
and concatenating explicitly fixes both failure modes.

`download.sh` runs standalone (sequential) on any machine; `download.slurm`
is an optional wrapper for Slurm clusters only — most users won't need it.
Before submitting a large array job, run one task end to end first and check
the result: correct byte count, `gzip -t` passes, real FASTQ content on the
first line.

## Repository layout

```
ribohunter/
  common.py            shared NCBI/EBI HTTP helpers
  discover.py           steps 1-2: search + per-sample keyword confirmation
  dedup.py               step 3: cross-listed sample deduplication
  verify_species.py      step 4: per-sample organism re-verification
  extract_metadata.py    step 5: cell line/tissue, treatment, adapter
  pair_rnaseq.py          step 6: RNA-seq pairing
  resolve_sra.py           step 7: SRX resolution + ENA run-level sizes
  build_manifest.py        assembles the multi-run-safe download manifest
  download.sh               downloader (sequential by default)
  download.slurm            optional Slurm array wrapper
examples/
  drosophila_melanogaster/  worked example: final table, manifest, contamination log
```

## Citation

If you use ribohunter in your work, please cite it:

```bibtex
@software{ribohunter,
  author = {Diaz de la Vega, Alejandro},
  title = {ribohunter: automated discovery, verification, and download of public ribosome profiling (Ribo-seq) datasets},
  year = {2026},
  url = {https://github.com/AlejandroDVG/Ribo_hunter}
}
```

A machine-readable [CITATION.cff](CITATION.cff) is also included, so GitHub's
"Cite this repository" button will generate this automatically.

## License

MIT — see [LICENSE](LICENSE).
