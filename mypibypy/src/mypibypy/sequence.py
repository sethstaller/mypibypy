"""Sequence utilities: fetching genomic DNA and basic biophysics.

Kept dependency-free (stdlib only) so the core is portable. Genome sequence is
pulled from the UCSC REST API (https://api.genome.ucsc.edu).
"""
from __future__ import annotations

import json
import math
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

UCSC_API = "https://api.genome.ucsc.edu"
ENSEMBL_REST = "https://rest.ensembl.org"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    """Reverse complement, preserving case."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def max_homopolymer(seq: str) -> int:
    """Length of the longest run of a single base."""
    seq = seq.upper()
    best = run = 0
    prev = ""
    for b in seq:
        run = run + 1 if b == prev else 1
        prev = b
        best = max(best, run)
    return best


def max_dinucleotide_repeat(seq: str) -> int:
    """Length (in bases) of the longest perfect dinucleotide repeat, e.g. ATATAT -> 6."""
    seq = seq.upper()
    n = len(seq)
    best = 0
    for start in range(n - 1):
        unit = seq[start : start + 2]
        if unit[0] == unit[1]:
            continue  # that's a homopolymer, handled elsewhere
        length = 2
        j = start + 2
        while j + 1 < n and seq[j : j + 2] == unit:
            length += 2
            j += 2
        best = max(best, length)
    return best


# --- Nearest-neighbor Tm (SantaLucia 1998, unified parameters) -------------
# ΔH kcal/mol, ΔS cal/(mol·K) for the 10 unique Watson-Crick NN pairs.
_NN_DH = {
    "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5,
    "GT": -8.4, "AC": -8.4, "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
    "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0,
}
_NN_DS = {
    "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7,
    "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9,
}
_R = 1.987  # cal/(mol·K)


def nn_tm(seq: str, na_molar: float = 0.33, oligo_molar: float = 0.25e-6) -> float:
    """Nearest-neighbor melting temperature (°C).

    Defaults: [Na+] = 0.33 M (~2x SSC), oligo 0.25 µM, no formamide. These are a
    fixed reference for *comparing* candidate windows, not a wash-condition model;
    add ~ -0.6 °C per % formamide for an in-buffer estimate.
    """
    s = seq.upper()
    if len(s) < 2 or any(b not in "ACGT" for b in s):
        return float("nan")
    dh = ds = 0.0
    # initiation (terminal bp)
    for end in (s[0], s[-1]):
        if end in "GC":
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1
    for i in range(len(s) - 1):
        pair = s[i : i + 2]
        dh += _NN_DH[pair]
        ds += _NN_DS[pair]
    # salt correction
    ds += 0.368 * (len(s) - 1) * math.log(na_molar)
    tm = (dh * 1000.0) / (ds + _R * math.log(oligo_molar / 4.0)) - 273.15
    return tm


@dataclass(frozen=True)
class GenomicInterval:
    build: str
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        return f"{self.build} {self.chrom}:{self.start:,}-{self.end:,}"


def fetch_sequence(interval: GenomicInterval) -> str:
    """Fetch the plus-strand sequence for a 1-based inclusive interval from UCSC."""
    # UCSC REST is 0-based half-open.
    url = (
        f"{UCSC_API}/getData/sequence?genome={interval.build}"
        f";chrom={interval.chrom};start={interval.start - 1};end={interval.end}"
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
    seq = data["dna"].upper()
    if len(seq) != interval.length:
        raise ValueError(
            f"UCSC returned {len(seq)} nt, expected {interval.length} for {interval}"
        )
    return seq


def rmsk_hits(interval: GenomicInterval) -> list[dict]:
    """RepeatMasker features overlapping the interval (empty list == repeat-free)."""
    url = (
        f"{UCSC_API}/getData/track?genome={interval.build};track=rmsk"
        f";chrom={interval.chrom};start={interval.start - 1};end={interval.end}"
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
    return data.get("rmsk", [])


def mask_repeats(seq: str, interval: GenomicInterval, hits: list[dict]) -> str:
    """Soft-mask (lowercase) positions overlapping RepeatMasker features.

    Window scanning rejects any window containing a non-ACGT (lowercase) base, so
    soft-masked repeats are excluded from probe placement.
    """
    chars = list(seq)
    offset = interval.start - 1  # genomic 0-based -> region index
    for h in hits:
        a = max(0, int(h["genoStart"]) - offset)
        b = min(len(chars), int(h["genoEnd"]) - offset)
        for i in range(a, b):
            chars[i] = chars[i].lower()
    return "".join(chars)


# --- Transcript sequence by accession --------------------------------------
#
# RNA piFISH targets a transcript, not an interval, so the input is an accession
# rather than coordinates. Two sources, because the lab's records use both:
# Ensembl/GENCODE ENST... and RefSeq NM_/NR_/XM_/XR_.
#
# Both fetchers validate what comes back. A truncated HTTP response that still
# parses as sequence is the quiet failure here: it would silently shorten a
# transcript, move every downstream offset, and design probes against the wrong
# place. `fetch_sequence` checks the returned length against the interval it
# asked for; these do the equivalent with what each API can tell us.

_ENSEMBL_PREFIXES = ("ENST", "ENSMUST", "ENSG", "ENSMUSG")
_REFSEQ_PREFIXES = ("NM_", "NR_", "XM_", "XR_", "NC_")


class TranscriptFetchError(RuntimeError):
    """A transcript could not be fetched, or what came back was not usable."""


#: HTTP statuses worth retrying. Ensembl REST and NCBI E-utilities both emit
#: transient 500s and rate-limit with 429 under normal use.
_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})


def _get(url: str, headers: dict | None = None, timeout: int = 30,
         retries: int = 4) -> bytes:
    """GET with retries on the transient statuses these services actually emit."""
    req = urllib.request.Request(url, headers={"User-Agent": "mypibypy", **(headers or {})})
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code in _RETRY_STATUSES and attempt < retries - 1:
                # Linear backoff. These are public shared services; retrying
                # harder than this is rude and would not help anyway.
                time.sleep(2 * (attempt + 1))
                continue
            raise TranscriptFetchError(
                f"{url} returned HTTP {exc.code}: {exc.reason}") from exc
        except urllib.error.URLError as exc:
            if isinstance(exc.reason, ssl.SSLCertVerificationError):
                # A TLS-inspecting proxy presents a root that is in the macOS
                # keychain but not in this interpreter's CA bundle, so curl
                # works and Python does not. Never paper over this by disabling
                # verification -- point OpenSSL at the system roots instead.
                raise TranscriptFetchError(
                    f"TLS certificate verification failed for {url}. This usually "
                    f"means a TLS-inspecting proxy whose root certificate is in "
                    f"the system keychain but not in this interpreter's CA "
                    f"bundle. Export the system roots and point OpenSSL at them:\n"
                    f"  security find-certificate -a -p /System/Library/Keychains"
                    f"/SystemRootCertificates.keychain > roots.pem\n"
                    f"  security find-certificate -a -p /Library/Keychains"
                    f"/System.keychain >> roots.pem\n"
                    f"  export SSL_CERT_FILE=$PWD/roots.pem\n"
                    f"Underlying error: {exc.reason}") from exc
            last = exc
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
                continue
            raise TranscriptFetchError(f"could not reach {url}: {exc.reason}") from exc
    raise TranscriptFetchError(f"{url} failed after {retries} attempts: {last}")


def _validate_sequence(seq: str, what: str) -> str:
    seq = seq.strip().upper().replace("U", "T")
    if not seq:
        raise TranscriptFetchError(f"{what}: empty sequence returned")
    bad = set(seq) - set("ACGTN")
    if bad:
        raise TranscriptFetchError(f"{what}: non-DNA characters {sorted(bad)[:5]}")
    return seq


def strip_version(accession: str) -> str:
    """``NM_000546.6`` -> ``NM_000546``. Ensembl lookups want the unversioned ID."""
    return accession.split(".", 1)[0]


def fetch_ensembl_cdna(accession: str, timeout: int = 30) -> str:
    """cDNA for an Ensembl transcript ID, via the Ensembl REST API.

    JSON is requested rather than FASTA on purpose: a truncated JSON body fails
    to parse, which turns a partial download into a loud error instead of a
    short sequence. Note the API wants an ``Accept`` header -- sending
    ``Content-Type`` instead returns 200 with an empty body, which is exactly
    the silent failure worth avoiding.
    """
    ident = strip_version(accession)
    url = f"{ENSEMBL_REST}/sequence/id/{ident}?type=cdna"
    raw = _get(url, headers={"Accept": "application/json"}, timeout=timeout)
    if not raw.strip():
        raise TranscriptFetchError(
            f"Ensembl returned an empty body for {accession} (url {url})")
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise TranscriptFetchError(
            f"Ensembl response for {accession} is not valid JSON, so it is "
            f"probably truncated: {exc}") from exc
    returned = str(data.get("id", ""))
    if returned and strip_version(returned) != ident:
        raise TranscriptFetchError(
            f"asked Ensembl for {ident} and got {returned} back")
    return _validate_sequence(str(data.get("seq", "")), f"Ensembl {accession}")


def refseq_length(accession: str, timeout: int = 30) -> int | None:
    """Recorded sequence length for a RefSeq accession, via esummary.

    Fetched separately so the downloaded sequence can be checked against an
    independently reported length. Returns None when esummary does not report
    one, rather than inventing a number to compare against.
    """
    url = (f"{NCBI_EUTILS}/esummary.fcgi?db=nuccore&id="
           f"{urllib.parse.quote(accession)}&retmode=json")
    try:
        data = json.loads(_get(url, timeout=timeout))
    except json.JSONDecodeError:
        return None
    result = data.get("result", {})
    for key, rec in result.items():
        if key == "uids" or not isinstance(rec, dict):
            continue
        slen = rec.get("slen")
        if isinstance(slen, int) and slen > 0:
            return slen
    return None


def fetch_refseq_transcript(accession: str, timeout: int = 30) -> str:
    """Sequence for a RefSeq accession, cross-checked against its recorded length.

    efetch hands back FASTA with no length in it, so truncation is invisible
    from the response alone. esummary reports the length independently, and the
    two are compared -- which is the closest equivalent here to what
    `fetch_sequence` does with an interval.
    """
    expected = refseq_length(accession, timeout=timeout)
    url = (f"{NCBI_EUTILS}/efetch.fcgi?db=nuccore&id="
           f"{urllib.parse.quote(accession)}&rettype=fasta&retmode=text")
    text = _get(url, timeout=timeout).decode("utf-8", "replace")
    lines = [ln.strip() for ln in text.splitlines()]
    if not lines or not lines[0].startswith(">"):
        raise TranscriptFetchError(
            f"NCBI did not return FASTA for {accession}: {text[:120]!r}")
    seq = _validate_sequence("".join(lines[1:]), f"RefSeq {accession}")
    if expected is not None and len(seq) != expected:
        raise TranscriptFetchError(
            f"RefSeq {accession}: downloaded {len(seq)} nt but esummary reports "
            f"{expected} nt -- the response is truncated or the accession is "
            f"ambiguous; refusing to design against it")
    return seq


def fetch_transcript(accession: str, timeout: int = 30) -> str:
    """Transcript sequence for an Ensembl or RefSeq accession.

    Dispatches on the accession prefix. Raises `TranscriptFetchError` for an
    unrecognised one rather than guessing a source, since guessing wrong returns
    a plausible sequence for the wrong molecule.
    """
    acc = accession.strip()
    upper = acc.upper()
    if upper.startswith(_ENSEMBL_PREFIXES):
        return fetch_ensembl_cdna(acc, timeout=timeout)
    if upper.startswith(_REFSEQ_PREFIXES):
        return fetch_refseq_transcript(acc, timeout=timeout)
    raise TranscriptFetchError(
        f"unrecognised accession {accession!r}: expected an Ensembl ID "
        f"({'/'.join(_ENSEMBL_PREFIXES)}) or a RefSeq accession "
        f"({'/'.join(_REFSEQ_PREFIXES)})")
