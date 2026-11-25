# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Remote lookups: Open Targets Platform, Ensembl (VEP, genes in region), STRING, g:Profiler.

All are open APIs without authentication. Each function returns plain JSON-able data
and is small enough for an agent's context.
"""

from __future__ import annotations

from typing import Any

from genor.errors import GenorError
from genor.lookup.http import request

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
ENSEMBL = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
STRING_URL = "https://version-12-0.string-db.org/api"
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


# ----------------------------------------------------------------- Open Targets


def ot_query(query: str, variables: dict[str, Any]) -> dict[str, Any]:
    out = request("POST", OT_URL, json={"query": query, "variables": variables})
    if out.get("errors"):
        raise GenorError(f"Open Targets: {out['errors'][0].get('message')}", code="remote_bad_request")
    return out["data"]


def ot_search(term: str, entity: str | None = None, size: int = 5) -> list[dict[str, Any]]:
    q = """query($q:String!,$e:[String!],$n:Int!){search(queryString:$q,entityNames:$e,page:{index:0,size:$n}){
      hits{id entity name description}}}"""
    data = ot_query(q, {"q": term, "e": [entity] if entity else None, "n": size})
    return data["search"]["hits"]


def ot_target(ensembl_id: str, n_diseases: int = 10) -> dict[str, Any]:
    q = """query($id:String!,$n:Int!){target(ensemblId:$id){id approvedSymbol approvedName biotype
      functionDescriptions tractability{modality label value}
      associatedDiseases(page:{index:0,size:$n}){count rows{score disease{id name}}}}}"""
    t = ot_query(q, {"id": ensembl_id, "n": n_diseases})["target"]
    if t is None:
        raise GenorError(f"No Open Targets target {ensembl_id}", code="not_found")
    t["tractability"] = [x for x in (t.get("tractability") or []) if x.get("value")]
    t["functionDescriptions"] = (t.get("functionDescriptions") or [])[:2]
    return t


def ot_disease_targets(efo_id: str, n: int = 25) -> dict[str, Any]:
    q = """query($id:String!,$n:Int!){disease(efoId:$id){id name associatedTargets(page:{index:0,size:$n}){
      count rows{score target{id approvedSymbol}}}}}"""
    d = ot_query(q, {"id": efo_id, "n": n})["disease"]
    if d is None:
        raise GenorError(f"No Open Targets disease {efo_id}", code="not_found")
    return d


def ot_variant_credible_sets(variant_id: str, n: int = 10) -> dict[str, Any]:
    """variant_id in Open Targets form CHR_POS_REF_ALT on GRCh38, e.g. 16_53767042_T_C."""
    q = """query($id:String!,$n:Int!){variant(variantId:$id){id rsIds
      credibleSets(page:{index:0,size:$n}){count rows{studyLocusId pValueMantissa pValueExponent
        study{id traitFromSource studyType}
        l2GPredictions(page:{index:0,size:3}){rows{score target{id approvedSymbol}}}}}}}"""
    v = ot_query(q, {"id": variant_id, "n": n})["variant"]
    if v is None:
        raise GenorError(f"Variant {variant_id} not in Open Targets", code="not_found",
                         hint="Use GRCh38 CHR_POS_REF_ALT. Open Targets is GRCh38-only.")
    v["note"] = (f"First {n} of {v['credibleSets']['count']} credible sets in API order (not ranked "
                 "by p-value); raise n to see more.")
    return v


def ot_study_locus_coloc(study_locus_id: str, n: int = 20) -> dict[str, Any]:
    q = """query($id:String!,$n:Int!){credibleSet(studyLocusId:$id){studyLocusId
      study{id traitFromSource}
      colocalisation(page:{index:0,size:$n}){count rows{colocalisationMethod h4 clpp
        otherStudyLocus{studyLocusId study{id studyType traitFromSource} }}}}}"""
    cs = ot_query(q, {"id": study_locus_id, "n": n})["credibleSet"]
    if cs is None:
        raise GenorError(f"No Open Targets credible set {study_locus_id}", code="not_found")
    return cs


# ----------------------------------------------------------------- Ensembl


def _ensembl(build: str) -> tuple[str, str]:
    from genor.guards import norm_build

    b = norm_build(build)
    if b is None:
        raise GenorError(f"Unknown genome build '{build}'", code="bad_param", hint="Use GRCh37 or GRCh38.")
    return b, ENSEMBL[b]


def vep(rsids: list[str], build: str = "GRCh38") -> list[dict[str, Any]]:
    if len(rsids) > 200:
        raise GenorError("VEP lookups are limited to 200 variants per call.", code="too_many")
    _, url = _ensembl(build)
    out = request("POST", f"{url}/vep/human/id",
                  json={"ids": rsids}, headers={"Content-Type": "application/json",
                                                "Accept": "application/json"})
    rows = []
    for v in out:
        tc = v.get("transcript_consequences") or []
        genes = sorted({t.get("gene_symbol") for t in tc if t.get("gene_symbol")})
        rows.append({"id": v.get("id"), "most_severe_consequence": v.get("most_severe_consequence"),
                     "location": f"{v.get('seq_region_name')}:{v.get('start')}",
                     "genes": genes[:10]})
    return rows


def genes_in_region(chrom: str, start: int, end: int, build: str = "GRCh37") -> list[dict[str, Any]]:
    if end - start > 5_000_000:
        raise GenorError("Region larger than 5 Mb.", code="region_too_large")
    build, url = _ensembl(build)
    try:
        out = request("GET", f"{url}/overlap/region/human/{chrom}:{start}-{end}",
                      params={"feature": "gene"}, headers={"Accept": "application/json"})
    except GenorError as e:
        if build == "GRCh37" and e.code == "remote_unavailable":
            e.hint = ("The GRCh37 Ensembl REST server (grch37.rest.ensembl.org) is unavailable. Lift the "
                      "region to GRCh38 and retry with --build GRCh38.")
        raise
    return [{"gene_id": g["gene_id"], "symbol": g.get("external_name"), "biotype": g.get("biotype"),
             "start": g["start"], "end": g["end"], "strand": g["strand"]} for g in out]


# ----------------------------------------------------------------- STRING


def string_network(genes: list[str], species: int = 9606, required_score: int = 700) -> dict[str, Any]:
    if len(genes) > 500:
        raise GenorError("STRING lookups are limited to 500 genes.", code="too_many")
    edges = request("POST", f"{STRING_URL}/json/network",
                    data={"identifiers": "\r".join(genes), "species": species,
                          "required_score": required_score, "caller_identity": "genor"})
    edges = [{"a": e["preferredName_A"], "b": e["preferredName_B"], "score": e["score"]} for e in edges]
    ppi = request("POST", f"{STRING_URL}/json/ppi_enrichment",
                  data={"identifiers": "\r".join(genes), "species": species,
                        "required_score": required_score, "caller_identity": "genor"})
    return {"string_version": "12.0", "n_edges": len(edges), "edges": edges[:300],
            "ppi_enrichment": ppi[0] if ppi else None}


# ----------------------------------------------------------------- g:Profiler


def gprofiler(genes: list[str], background: list[str] | None = None,
              sources: list[str] | None = None, threshold: float = 0.05) -> dict[str, Any]:
    body: dict[str, Any] = {"organism": "hsapiens", "query": genes, "user_threshold": threshold,
                            "significance_threshold_method": "g_SCS", "no_evidences": True,
                            "sources": sources or ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "HPA"]}
    if background:
        body["domain_scope"] = "custom"
        body["background"] = background
    out = request("POST", GPROFILER_URL, json=body)
    res = [{"source": r["source"], "term_id": r["native"], "name": r["name"],
            "p_value": r["p_value"], "term_size": r["term_size"],
            "intersection_size": r["intersection_size"]} for r in out.get("result", [])]
    meta = out.get("meta", {})
    return {"n_terms": len(res), "terms": res[:100],
            "gprofiler_version": meta.get("version"),
            "background": "custom" if background else "annotated genes (default)"}
