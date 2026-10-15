"""Genome-wide specificity check.

UCSC web BLAT is behind a Cloudflare bot-check, so this falls back to NCBI's
URL BLAST API (programmatic, no captcha, but slow: submit -> poll -> parse).
Use against the mouse genome to confirm each split-probe footprint and its 25-nt
pads are single-copy.
"""
from __future__ import annotations

import time
import urllib.error
import urllib.parse
import urllib.request

NCBI = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _post(params: dict, retries: int = 4) -> str:
    data = urllib.parse.urlencode(params).encode()
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                NCBI, data=data, headers={"User-Agent": "mypibypy/0.1"})
            return urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:  # 502/503 from NCBI are common + transient
            last = e
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(10 * (attempt + 1))
                continue
            raise
    raise last


def submit(fasta: str, organism_taxid: str = "10090") -> tuple[str, int]:
    """Submit a (multi-)FASTA blastn job against mouse genomic sequence; return (RID, RTOE)."""
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "refseq_genomes",
        "ENTREZ_QUERY": f"txid{organism_taxid}[ORGN]",
        "QUERY": fasta,
        "WORD_SIZE": "11",          # sensitive (short queries)
        "HITLIST_SIZE": "20",
    }
    html = _post(params)
    rid = rtoe = None
    for line in html.splitlines():
        line = line.strip()
        if line.startswith("RID = "):
            rid = line[6:].strip()
        elif line.startswith("RTOE = "):
            rtoe = int(line[7:].strip())
    if not rid:
        raise RuntimeError("NCBI BLAST did not return an RID")
    return rid, (rtoe or 30)


def poll(rid: str, interval: int = 20, timeout: int = 900) -> str:
    """Poll until the job is READY; return Text report."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        html = _post({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"})
        if "Status=READY" in html:
            return _post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text"})
        if "Status=FAILED" in html or "Status=UNKNOWN" in html:
            raise RuntimeError(f"NCBI BLAST job {rid} failed")
        time.sleep(interval)
    raise TimeoutError(f"NCBI BLAST job {rid} timed out")


def to_fasta(named_seqs: dict[str, str]) -> str:
    return "\n".join(f">{name}\n{seq}" for name, seq in named_seqs.items())
