"""Cross-hybridization against the amplifier pool and the rest of the order.

General probe-design tools score each probe in isolation: genome hits, folding,
self-complementarity. None of them knows that a piFISH primary is ordered into a
*pool* that already contains P2/P3/P4 for every channel. This module asks the
questions that only make sense in that context.

The piFISH coincidence mechanism, verified against the validated systems in
`systems.py` rather than assumed:

    P1L = arm_left + rc(left25)          P1R = rc(right25) + arm_right

    arm_left  = [ top_L: 14 nt ][ mid_L: 8 nt ]
    arm_right = [ mid_R: 8 nt ][ top_R: 14 nt ]

    P2 anneal region (29 nt) = rc(top_R) + 1 nt + rc(top_L)

So P2's landing site exists only when both primaries are bound to adjacent
target: each contributes 14 nt of it. mid_L and mid_R form a short (4 bp) stem
across the junction -- the "pi bond" that holds the pair together.

That geometry defines the failure we care about, which single-probe QC cannot
see: **if a primary's target pad happens to extend its own 14-nt half-site, one
primary alone can recruit P2** and the coincidence requirement is lost. That is
`lone_primary_nucleation`, and it is the check worth having.

The other checks (pad vs amplifier repeats, primary vs the other channel's arms,
primary vs primary within an order) are cheaper and more routine, but still
invisible to per-probe tools.

Thresholds are calibrated against the three existing validated designs -- see
`scripts/calibrate_crosshyb.py`. They are starting points, not measured limits.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .sequence import _NN_DH, _NN_DS, revcomp
from .systems import AmplifierSystem

#: Geometry of an amplifier arm, from the paper and confirmed in `systems.py`.
TOP_LEN = 14
MID_LEN = 8

#: A perfect complementary run at or above this length is reported.
REPORT_RUN = 8

#: Default flag thresholds (nt of contiguous complementarity).
#: Calibrated so the three validated designs pass -- see module docstring.
MAX_PAD_VS_AMPLIFIER = 11
MAX_PRIMARY_VS_ARMS = 12
MAX_PRIMARY_VS_PRIMARY = 12

#: A lone primary presents TOP_LEN nt of the P2 site by design, and the 8-nt
#: middle region sits between that top and the pad -- so a *contiguous* duplex
#: can never run past the top. The real lone-primary risk is the pad pairing
#: with the P2 anneal region as a SECOND segment: top (14) + pad (n) together
#: can hold P2 without the partner primary. This caps that second segment.
MAX_NUCLEATION_PAD_CONTRIBUTION = 7


def _clean(seq: str) -> str:
    """Strip IDT modification tags like /5Cy5/ and uppercase."""
    out, i = [], 0
    while i < len(seq):
        if seq[i] == "/":
            j = seq.find("/", i + 1)
            if j == -1:
                break
            i = j + 1
            continue
        out.append(seq[i])
        i += 1
    return "".join(out).upper()


def duplex_dg(seq: str, temp_c: float = 40.0) -> float:
    """Gibbs free energy (kcal/mol) of a perfect duplex of `seq` at `temp_c`.

    Nearest-neighbour dH/dS at 1 M Na+ with helix initiation. Negative is stable.
    Salt and formamide are not applied here: this is used to *rank* candidate
    spurious duplexes, and `thermo` handles buffer-accurate numbers.
    """
    s = seq.upper()
    if len(s) < 2 or any(b not in "ACGT" for b in s):
        return 0.0
    dh = ds = 0.0
    for end in (s[0], s[-1]):
        if end in "GC":
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1
    for i in range(len(s) - 1):
        dh += _NN_DH[s[i : i + 2]]
        ds += _NN_DS[s[i : i + 2]]
    return dh - (temp_c + 273.15) * ds / 1000.0


@dataclass(frozen=True)
class Duplex:
    """A gapless stretch where `a` can pair with `b`."""

    length: int
    dg: float
    a_start: int          # 0-based start of the run in a
    b_start: int          # 0-based start in b (on b's own 5'->3' coordinates)
    matched: str          # the run as it reads on a, 5'->3'

    def __bool__(self) -> bool:
        return self.length > 0


EMPTY = Duplex(0, 0.0, -1, -1, "")


def best_duplex(a: str, b: str, temp_c: float = 40.0) -> Duplex:
    """Longest perfect antiparallel complementary run between `a` and `b`.

    Equivalent to the longest common substring of `a` and revcomp(`b`). Ties are
    broken by stability (most negative dG), so a GC-rich run outranks an
    equally long AT-rich one. Mismatched duplexes are handled implicitly: a
    mismatch terminates a run, and the surviving runs are what can nucleate.
    """
    a, b = _clean(a), _clean(b)
    rb = revcomp(b)
    if not a or not rb:
        return EMPTY
    best = EMPTY
    prev = [0] * (len(rb) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(rb) + 1)
        ai = a[i - 1]
        for j in range(1, len(rb) + 1):
            if ai != rb[j - 1]:
                continue
            n = prev[j - 1] + 1
            cur[j] = n
            if n >= best.length and n >= 2:
                run = a[i - n : i]
                dg = duplex_dg(run, temp_c)
                if n > best.length or dg < best.dg:
                    # rb index j-n maps back to b: rb[k] == revcomp(b)[k]
                    best = Duplex(n, dg, i - n, len(b) - (j - n) - n, run)
        prev = cur
    return best


# --- Amplifier pool decomposition -----------------------------------------

@dataclass(frozen=True)
class PoolPart:
    """One functionally distinct piece of the amplifier pool."""

    name: str
    seq: str
    role: str   # "anneal" | "repeat" | "signal" | "arm"


def arm_regions(system: AmplifierSystem) -> dict[str, str]:
    """Split both arms into their top (P2-binding) and middle (pi-bond) regions."""
    return {
        "top_L": system.arm_left[:TOP_LEN],
        "mid_L": system.arm_left[TOP_LEN:],
        "mid_R": system.arm_right[:MID_LEN],
        "top_R": system.arm_right[MID_LEN:],
    }


def p2_anneal(system: AmplifierSystem) -> tuple[str, int]:
    """Locate P2's central anneal region. Returns (sequence, offset in P2).

    Found by searching for the reverse complements of the two arm tops rather
    than by hardcoding an offset, so this stays correct for new systems.
    Raises if the expected architecture is not present.
    """
    r = arm_regions(system)
    want_r, want_l = revcomp(r["top_R"]), revcomp(r["top_L"])
    i = system.p2.find(want_r)
    j = system.p2.find(want_l)
    if i == -1 or j == -1:
        raise ValueError(
            f"system {system.name}: P2 does not contain rc(top_R) and rc(top_L); "
            "the arm/P2 architecture is not what this module assumes"
        )
    if j < i:
        i, j = j, i
    return system.p2[i : j + TOP_LEN], i


def repeat_unit(seq: str, min_len: int = 10, max_len: int = 40) -> str | None:
    """Smallest tandem repeat unit that tiles the start of `seq`, if any."""
    s = _clean(seq)
    for n in range(min_len, min(max_len, len(s) // 2) + 1):
        unit = s[:n]
        if s[n : 2 * n] == unit:
            return unit
    return None


def pool_parts(systems: list[AmplifierSystem]) -> list[PoolPart]:
    """Every functionally distinct sequence a new primary will share a tube with."""
    parts: list[PoolPart] = []
    for s in systems:
        try:
            anneal, _ = p2_anneal(s)
            parts.append(PoolPart(f"{s.name}:P2_anneal", anneal, "anneal"))
        except ValueError:
            parts.append(PoolPart(f"{s.name}:P2_whole", s.p2, "anneal"))
        for label, seq in (("P2", s.p2), ("P3", s.p3)):
            unit = repeat_unit(seq)
            parts.append(PoolPart(f"{s.name}:{label}_repeat", unit or seq,
                                  "repeat" if unit else "anneal"))
        parts.append(PoolPart(s.p4_name, _clean(s.p4), "signal"))
        for name, seq in arm_regions(s).items():
            parts.append(PoolPart(f"{s.name}:{name}", seq, "arm"))
    return parts


# --- Checks ----------------------------------------------------------------

@dataclass
class Finding:
    kind: str
    subject: str
    against: str
    duplex: Duplex
    threshold: int
    note: str = ""

    @property
    def flagged(self) -> bool:
        return self.duplex.length > self.threshold

    def __str__(self) -> str:
        mark = "FLAG" if self.flagged else "ok  "
        return (f"{mark} {self.kind:<26} {self.subject:<22} vs {self.against:<26} "
                f"{self.duplex.length:>2} bp  dG {self.duplex.dg:>6.1f}  "
                f"{self.duplex.matched}{('  ' + self.note) if self.note else ''}")


@dataclass
class CrossHybReport:
    findings: list[Finding] = field(default_factory=list)

    @property
    def flagged(self) -> list[Finding]:
        return [f for f in self.findings if f.flagged]

    @property
    def ok(self) -> bool:
        return not self.flagged

    def worst(self, kind: str | None = None) -> Finding | None:
        pool = [f for f in self.findings if kind is None or f.kind == kind]
        return max(pool, key=lambda f: f.duplex.length, default=None)

    def summary(self) -> str:
        if self.ok:
            w = self.worst()
            return (f"no cross-hybridization flags "
                    f"(worst {w.duplex.length} bp, {w.kind})" if w else "no checks run")
        return f"{len(self.flagged)} cross-hybridization flag(s): " + "; ".join(
            f"{f.subject} vs {f.against} {f.duplex.length} bp" for f in self.flagged
        )


def _pad_of(primary: str, system: AmplifierSystem) -> str:
    """The target-binding pad of a primary, whichever side its arm is on."""
    s = primary.upper()
    if s.startswith(system.arm_left):
        return s[len(system.arm_left) :]
    if s.endswith(system.arm_right):
        return s[: -len(system.arm_right)]
    return s


def lone_primary_nucleation(primary: str, primary_name: str, system: AmplifierSystem,
                            side: str, temp_c: float = 40.0) -> Finding:
    """Can this primary alone hold P2, without its partner?

    This is the check that matters, and the one no per-probe tool can make.

    By design a lone primary offers exactly TOP_LEN nt of P2's 29-nt anneal
    region. The 8-nt middle region lies between that top and the target pad, and
    is deliberately not complementary to P2 -- so a lone primary's duplex cannot
    simply run on past the top.

    The way it *can* fail is the pad pairing with the anneal region as a second
    segment. Top (14 bp) plus a pad contribution of n bp, bridged by a bulge, can
    be enough to recruit P2 with only one primary bound -- which is exactly the
    coincidence requirement that gives piFISH its specificity. So the flag is on
    the pad's own contribution to the anneal region, with the whole-primary
    duplex reported alongside for context.
    """
    try:
        anneal, _ = p2_anneal(system)
    except ValueError as exc:
        return Finding("lone_primary_nucleation", f"{primary_name} ({side})",
                       f"{system.name}:P2", EMPTY,
                       MAX_NUCLEATION_PAD_CONTRIBUTION, f"not checked: {exc}")
    pad = _pad_of(primary, system)
    whole = best_duplex(primary, anneal, temp_c)
    pad_d = best_duplex(pad, anneal, temp_c)
    note = (f"top contributes {whole.length} bp (by design {TOP_LEN}); "
            f"pad adds {pad_d.length} bp -> {whole.length + pad_d.length} bp total "
            f"available to a lone primary")
    return Finding("lone_primary_nucleation", f"{primary_name} ({side})",
                   f"{system.name}:P2_anneal", pad_d,
                   MAX_NUCLEATION_PAD_CONTRIBUTION, note)


def check_pair(p1l: str, p1r: str, system: AmplifierSystem,
               systems: list[AmplifierSystem] | None = None,
               p1l_name: str = "P1L", p1r_name: str = "P1R",
               temp_c: float = 40.0) -> CrossHybReport:
    """Screen one primary pair against the whole amplifier pool.

    `systems` defaults to every registered system, because a multiplexed
    experiment puts both channels' amplifiers in the same sample.
    """
    from .systems import SYSTEMS

    systems = systems or list(dict.fromkeys(SYSTEMS.values()))
    rep = CrossHybReport()

    # The pads as synthesized: P1L carries its pad 3' of the arm, P1R 5' of it.
    pad_l = p1l[len(system.arm_left) :]
    pad_r = p1r[: -len(system.arm_right)] if system.arm_right else p1r

    # 1. Target pads against every amplifier part.
    for pad_name, pad in ((f"{p1l_name}.pad", pad_l), (f"{p1r_name}.pad", pad_r)):
        for part in pool_parts(systems):
            if part.role == "arm" and part.name.startswith(system.name):
                continue  # own arms are checked as nucleation, below
            d = best_duplex(pad, part.seq, temp_c)
            if d.length >= REPORT_RUN:
                rep.findings.append(
                    Finding("pad_vs_amplifier", pad_name, part.name, d,
                            MAX_PAD_VS_AMPLIFIER))

    # 2. Whole primary against the other channels' arms.
    for other in systems:
        if other.name == system.name:
            continue
        for arm_name, arm in (("arm_left", other.arm_left), ("arm_right", other.arm_right)):
            for pname, p in ((p1l_name, p1l), (p1r_name, p1r)):
                d = best_duplex(p, arm, temp_c)
                if d.length >= REPORT_RUN:
                    rep.findings.append(
                        Finding("primary_vs_other_arms", pname,
                                f"{other.name}:{arm_name}", d, MAX_PRIMARY_VS_ARMS))

    # 3. The two primaries against each other, beyond the designed pi bond.
    d = best_duplex(p1l, p1r, temp_c)
    rep.findings.append(
        Finding("primary_vs_partner", p1l_name, p1r_name, d, MAX_PRIMARY_VS_PRIMARY,
                "4 bp pi bond (GGTC/GACC) is by design"))

    # 4. Lone-primary P2 nucleation -- the coincidence-detection check.
    rep.findings.append(lone_primary_nucleation(p1l, p1l_name, system, "left", temp_c))
    rep.findings.append(lone_primary_nucleation(p1r, p1r_name, system, "right", temp_c))

    return rep


#: A pad duplex this close to full length means the two primaries target
#: opposite strands of the same site -- a sense/antisense pair, not an accident.
SENSE_ANTISENSE_FRACTION = 0.9


def check_order(primaries: dict[str, tuple[str, AmplifierSystem]],
                temp_c: float = 40.0) -> CrossHybReport:
    """Every primary in one order against every other.

    `primaries` maps name -> (sequence, its system). Pairs that belong to the
    same probe pair will show their designed pi bond; everything else should be
    well below threshold.

    Sense/antisense primaries for one target are reported separately: their pads
    are reverse complements by construction, so they pair over their full length.
    That is unavoidable, not a design error -- but it does mean the two compete
    in a shared tube, so it is surfaced rather than hidden.
    """
    rep = CrossHybReport()
    names = list(primaries)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            seq_a, sys_a = primaries[a]
            seq_b, sys_b = primaries[b]
            d = best_duplex(seq_a, seq_b, temp_c)
            if d.length < REPORT_RUN:
                continue
            pad_a = _pad_of(seq_a, sys_a)
            pad_b = _pad_of(seq_b, sys_b)
            shortest_pad = min(len(pad_a), len(pad_b)) or 1
            if d.length >= SENSE_ANTISENSE_FRACTION * shortest_pad:
                rep.findings.append(
                    Finding("sense_antisense_pair", a, b, d,
                            threshold=10 ** 6,  # informational: never a hard flag
                            note=("pads are reverse complements -- opposite strands of "
                                  "one target; they will dimerize if pooled together")))
            else:
                rep.findings.append(
                    Finding("order_cross_pairing", a, b, d, MAX_PRIMARY_VS_PRIMARY))
    for name, (seq, _) in primaries.items():
        d = best_duplex(seq, seq, temp_c)
        if d.length >= REPORT_RUN:
            rep.findings.append(
                Finding("self_complementarity", name, name, d, MAX_PRIMARY_VS_PRIMARY,
                        "self-dimer / hairpin potential"))
    return rep
