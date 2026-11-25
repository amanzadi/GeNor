# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Guard rails enforced in code.

- Summary statistics only: individual-level genotype inputs are refused.
- Genome build and ancestry must match the LD reference a method uses.
"""

from __future__ import annotations

import bz2
import gzip
import lzma
from pathlib import Path
from typing import TextIO

from genor.errors import GenorError, IndividualLevelDataError, NeedsUserDecision

POLICY = (
    "GeNor works on GWAS summary statistics only. It does not read individual-level "
    "genotype or phenotype data (PLINK .bed/.pgen/.bgen, VCFs with sample genotypes). "
    "Run your GWAS elsewhere and bring the summary statistics. Public LD reference panels "
    "(1000 Genomes) are handled internally by GeNor."
)

INDIVIDUAL_EXT = {".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam", ".bgen", ".sample",
                  ".ped", ".map", ".tped", ".tfam", ".gen", ".dosage", ".cram", ".bam"}

MHC = {"GRCh37": (6, 25_000_000, 34_000_000), "GRCh38": (6, 25_000_000, 34_000_000)}


def _suffixes(path: Path) -> list[str]:
    return [s.lower() for s in path.suffixes]


def refuse_individual_level(path: Path) -> None:
    sfx = [s for s in _suffixes(path) if s not in {".gz", ".bgz", ".zst", ".bz2", ".xz"}]
    if sfx and sfx[-1] in INDIVIDUAL_EXT:
        raise IndividualLevelDataError(
            f"{path.name} looks like individual-level genotype data ({sfx[-1]}).",
            code="individual_level_refused",
            hint=POLICY,
        )
    if any(s == ".vcf" for s in sfx) or path.name.lower().endswith((".vcf.gz", ".vcf.bgz", ".vcf")):
        _check_vcf(path)


def open_text(path: Path) -> TextIO:
    """Open a plain, gzip, bzip2 or xz text file; refuse zstd (no stdlib reader)."""
    name = path.name.lower()
    if name.endswith(".zst"):
        raise GenorError(f"{path.name} is zstd-compressed, which GeNor cannot read.",
                         code="unsupported_compression",
                         hint=f"Decompress it first (zstd -d {path.name}) or recompress with gzip.")
    kw = {"encoding": "utf-8", "errors": "replace"}
    if name.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", **kw)
    if name.endswith(".bz2"):
        return bz2.open(path, "rt", **kw)
    if name.endswith(".xz"):
        return lzma.open(path, "rt", **kw)
    return open(path, **kw)


def _check_vcf(path: Path) -> None:
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= 9:
                    return  # sites-only VCF
                first = next(f, "").split("\t")
                fmt = first[8].split(":") if len(first) > 8 else []
                if "ES" in fmt and ("SE" in fmt or "LP" in fmt):
                    if len(cols) > 10:
                        raise IndividualLevelDataError(
                            f"{path.name} is a GWAS-VCF with {len(cols) - 9} trait columns; "
                            "GeNor reads one trait per file.",
                            code="gwas_vcf_multi",
                            hint="Split the file (bcftools view -s <trait>) and load each trait.",
                        )
                    return  # GWAS-VCF (summary statistics): allowed
                if "GT" in fmt or "DS" in fmt or "GP" in fmt:
                    raise IndividualLevelDataError(
                        f"{path.name} contains per-sample genotypes ({len(cols) - 9} samples).",
                        code="individual_level_refused",
                        hint=POLICY,
                    )
                return
    return


def norm_build(build: str | int | None) -> str | None:
    if build is None:
        return None
    b = str(build).strip().lower().replace("grch", "").replace("hg", "").replace("b", "")
    out = {"37": "GRCh37", "19": "GRCh37", "38": "GRCh38"}.get(b)
    if out is None:
        raise GenorError(f"Unrecognised genome build '{build}'.", code="bad_build",
                         hint="Use GRCh37 (also hg19, b37) or GRCh38 (also hg38, b38). Older builds "
                              "(hg18/GRCh36) must be lifted to one of these outside GeNor.")
    return out


ANCESTRIES = {"EUR", "EAS", "AFR", "AMR", "SAS", "MID", "CSA", "multi"}


def norm_ancestry(a: str | None) -> str | None:
    if a is None:
        return None
    up = a.strip().upper()
    return "multi" if up in {"MULTI", "TRANS", "MIXED", "MULTI-ANCESTRY"} else up


def check_build(sumstats_build: str | None, panel_build: str, what: str) -> None:
    if sumstats_build is None:
        raise NeedsUserDecision(
            f"The genome build of the summary statistics is unknown, and {what} uses {panel_build}.",
            question="Which genome build are the summary statistics on (GRCh37 or GRCh38)?",
            options=["GRCh37", "GRCh38"],
        )
    if sumstats_build != panel_build:
        raise GenorError(
            f"Build mismatch: sumstats are {sumstats_build} but {what} is {panel_build}.",
            code="build_mismatch",
            hint=f"Lift the sumstats over first: genor sumstats liftover <run> --to {panel_build}",
        )


def check_ancestry(sumstats_anc: str | None, panel_anc: str, what: str) -> None:
    if sumstats_anc in (None, "UNKNOWN"):
        raise NeedsUserDecision(
            f"The ancestry of the GWAS sample is unknown, and {what} is {panel_anc}.",
            question="What is the (majority) genetic ancestry of the GWAS sample?",
            options=sorted(ANCESTRIES - {"multi"}),
        )
    if sumstats_anc == "multi":
        raise NeedsUserDecision(
            f"The GWAS is multi-ancestry; the {what} is single-ancestry ({panel_anc}).",
            question=(f"Use the {panel_anc} panel anyway (approximate LD; common practice when "
                      f"most of the sample is {panel_anc}), or stop?"),
            options=[f"use {panel_anc}", "stop"],
        )
    if sumstats_anc != panel_anc:
        raise GenorError(
            f"Ancestry mismatch: GWAS is {sumstats_anc} but {what} is {panel_anc}.",
            code="ancestry_mismatch",
            hint=f"Use a {sumstats_anc} reference (genor refdata list) or confirm with the user.",
        )
