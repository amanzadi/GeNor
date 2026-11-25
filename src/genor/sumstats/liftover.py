# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Native UCSC chain-file liftover for point positions (SNPs).

Positions are 1-based in and out. Variants that fall outside aligned blocks, map
to a different chromosome, or map to the reverse strand are dropped (and counted).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class _Blocks:
    t_start: np.ndarray  # 0-based source start of each aligned block
    t_end: np.ndarray
    q_chrom: np.ndarray
    q_start: np.ndarray  # 0-based target start
    q_strand: np.ndarray  # +1 / -1
    q_size: np.ndarray


# GeNor codes sex/mito chromosomes numerically (load.CHROM_MAP); chain files use names
TO_CHAIN = {"23": "X", "24": "Y", "26": "M"}
FROM_CHAIN = {"X": "23", "Y": "24", "M": "26", "MT": "26"}


def _clean(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def load_chain(path: Path) -> dict[str, _Blocks]:
    rows: dict[str, list[tuple[int, int, str, int, int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        t_chrom = q_chrom = ""
        t_pos = q_pos = 0
        q_strand = 1
        q_size = 0
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "chain":
                t_chrom = _clean(parts[2])
                t_pos = int(parts[5])
                q_chrom = _clean(parts[7])
                q_size = int(parts[8])
                q_strand = 1 if parts[9] == "+" else -1
                q_pos = int(parts[10])
                continue
            size = int(parts[0])
            rows.setdefault(t_chrom, []).append((t_pos, t_pos + size, q_chrom, q_pos, q_strand, q_size))
            if len(parts) == 3:
                t_pos += size + int(parts[1])
                q_pos += size + int(parts[2])
    out = {}
    for chrom, r in rows.items():
        r.sort(key=lambda x: x[0])
        a = list(zip(*r, strict=True))
        out[chrom] = _Blocks(np.array(a[0], dtype=np.int64), np.array(a[1], dtype=np.int64),
                             np.array(a[2], dtype=object), np.array(a[3], dtype=np.int64),
                             np.array(a[4], dtype=np.int8), np.array(a[5], dtype=np.int64))
    return out


def lift(chain: dict[str, _Blocks], chrom: pd.Series, pos: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Return (new_chrom, new_pos); unmapped entries are NA."""
    new_chrom = pd.Series(pd.NA, index=chrom.index, dtype="object")
    new_pos = pd.Series(pd.NA, index=chrom.index, dtype="Int64")
    chrom_s = chrom.astype(str)
    for c in chrom_s.unique():
        cc = TO_CHAIN.get(c, c)
        blocks = chain.get(cc)
        if blocks is None and cc == "M":
            cc, blocks = "MT", chain.get("MT")
        if blocks is None:
            continue
        idx = chrom_s.index[chrom_s == c]
        p0 = pos.loc[idx].to_numpy(dtype=np.int64) - 1
        k = np.searchsorted(blocks.t_start, p0, side="right") - 1
        ok = k >= 0
        kk = np.where(ok, k, 0)
        ok &= p0 < blocks.t_end[kk]
        off = p0 - blocks.t_start[kk]
        fwd = blocks.q_strand[kk] == 1
        q = np.where(fwd, blocks.q_start[kk] + off, blocks.q_size[kk] - (blocks.q_start[kk] + off) - 1)
        ok &= fwd  # reverse-strand mappings would need allele complementing; drop them
        ok &= blocks.q_chrom[kk] == cc
        new_chrom.loc[idx[ok]] = FROM_CHAIN.get(cc, cc)
        new_pos.loc[idx[ok]] = q[ok] + 1
    return new_chrom, new_pos
