"""Presentation helpers shared by the Streamlit app.

These live here rather than in `app.py` so they can be imported and tested: a
Streamlit script executes its page body on import, so anything defined inside it
is awkward to reach from a test. The app keeps the layout; this keeps the logic
that turns a `DesignResult` into rows and bytes.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from .workbook import write_workbook


def candidate_rows(result) -> list[dict]:
    """One row per ordered pair: composition, thermodynamics, and both screens.

    `accessibility` and `specificity` render the *absence* of a result as words
    rather than as a blank or a zero -- "not evaluated" and "NOT SCREENED" --
    because a table cell is exactly where a missing screen would otherwise look
    like a passing one.
    """
    rows = []
    for i, (pair, screen, access) in enumerate(result.rows(), start=1):
        cand = pair.candidate
        rows.append({
            "#": i,
            "footprint": str(cand.interval),
            "GC L": round(cand.left.gc, 3),
            "GC R": round(cand.right.gc, 3),
            "Tm L": round(cand.left.tm, 1),
            "Tm R": round(cand.right.tm, 1),
            "accessibility": (round(access.pair_nucleation_p, 3)
                              if access.evaluated else "not evaluated"),
            "specificity": screen.verdict if screen.screened else "NOT SCREENED",
            "off-target genes": ", ".join(screen.offtarget_genes[:4]) or "-",
            "QC flags": ", ".join(cand.flags) or "none",
            "P1L": pair.p1l_name,
            "P1R": pair.p1r_name,
        })
    return rows


def workbook_bytes(results, note: str) -> bytes:
    """Render the order workbook into memory.

    Via a temporary directory, never into the repository: the workbook carries
    probe sequences, and the browser download is the only place they should go.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "order.xlsx"
        write_workbook(str(path), results, note)
        return path.read_bytes()


def unscreened_pairs(result) -> int:
    """How many ordered pairs have no specificity result. Drives the UI warning."""
    return sum(1 for screen in result.pair_screens if not screen.screened)
