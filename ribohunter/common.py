"""Shared HTTP/NCBI/EBI helpers used across the ribohunter pipeline steps."""
import json
import time
import urllib.parse
import urllib.request

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# NCBI allows ~3 unauthenticated requests/second; stay comfortably under that.
NCBI_DELAY_S = 0.34


def http_get(url, timeout=30, retries=3, retry_delay=1.5):
    """GET a URL and return the decoded response body, retrying on transient errors."""
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ribohunter/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001 - broad on purpose, this is a best-effort retry loop
            last_err = e
            time.sleep(retry_delay * (attempt + 1))
    raise RuntimeError(f"GET failed after {retries} attempts: {url}") from last_err


def esearch(term, db="gds", retmax=100000):
    """Run an NCBI esearch query and return the list of matching UIDs."""
    params = {"db": db, "term": term, "retmax": str(retmax), "retmode": "json"}
    url = f"{NCBI_EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    data = json.loads(http_get(url))
    return data.get("esearchresult", {}).get("idlist", [])


def esearch_single_uid(accession, db="gds"):
    """Resolve a single GEO sample accession (a GSM) to its UID, or None if not found."""
    ids = esearch(f"{accession}[ACCN] AND GSM[Filter]", db=db)
    return ids[0] if ids else None


def esummary(uids, db="gds"):
    """Batch-fetch esummary records for a list of UIDs (chunked at 100 per request)."""
    result = {}
    for chunk in _chunks(uids, 100):
        params = {"db": db, "id": ",".join(chunk), "retmode": "json"}
        url = f"{NCBI_EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params)
        data = json.loads(http_get(url))
        chunk_result = data.get("result", {})
        chunk_result.pop("uids", None)  # NCBI's esummary result dict has a stray "uids" list key alongside the per-UID entries
        result.update(chunk_result)
        time.sleep(NCBI_DELAY_S)
    return result


def fetch_geo_full_record(accession, view="full"):
    """Fetch a GEO accession's plain-text record (GSM or GSE) and parse it into
    a dict of field -> list-of-values (GEO repeats keys like Sample_characteristics_ch1)."""
    params = {"acc": accession, "targ": "self", "form": "text", "view": view}
    url = GEO_ACC_CGI + "?" + urllib.parse.urlencode(params)
    text = http_get(url)
    fields = {}
    for line in text.splitlines():
        if not line.startswith("!"):
            continue
        line = line[1:]
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        fields.setdefault(key.strip(), []).append(val.strip())
    return fields


def ena_filereport(srx_accession, retries=3):
    """Resolve an SRA experiment accession to its run-level file report via ENA."""
    params = {
        "accession": srx_accession,
        "result": "read_run",
        "fields": "run_accession,fastq_bytes,fastq_ftp,sample_accession",
        "format": "json",
    }
    url = ENA_FILEREPORT + "?" + urllib.parse.urlencode(params)
    txt = http_get(url, retries=retries)
    if not txt.strip():
        return []
    return json.loads(txt)


def _chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]
