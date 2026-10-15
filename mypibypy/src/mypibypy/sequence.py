"""Sequence utilities: fetching genomic DNA and basic biophysics.

Kept dependency-free (stdlib only) so the core is portable. Genome sequence is
pulled from the UCSC REST API (https://api.genome.ucsc.edu).
"""
from __future__ import annotations

import json
import math
import urllib.request
from dataclasses import dataclass

UCSC_API = "https://api.genome.ucsc.edu"

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
