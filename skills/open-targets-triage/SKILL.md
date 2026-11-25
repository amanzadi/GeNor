---
name: open-targets-triage
description: Query the Open Targets Platform (GraphQL API; the former Genetics portal is merged into it) with genor lookup search|target|disease|variant|coloc to find public credible sets, locus-to-gene (L2G) scores, colocalisations, target-disease association scores and tractability, and to triage candidate drug targets. Use when the user asks about Open Targets, L2G, druggability or tractability, known disease-target associations, or wants to rank candidate genes from GWAS loci as targets.
---
# Open Targets triage

GeNor queries the Open Targets Platform GraphQL API (`api.platform.opentargets.org/api/v4/graphql`).
The old Open Targets Genetics portal has been retired and its content (GWAS credible sets, L2G,
colocalisation) now lives in the Platform. Everything in the Platform is **GRCh38**.

| Command | MCP | Returns |
|---|---|---|
| `genor lookup search TERM [--entity target\|disease\|variant\|study]` | `lookup_open_targets(kind="search", query=..., n=...)` | Hits with `id`, `entity`, `name`, `description` |
| `genor lookup target ENSEMBL_ID [--n 10]` | `kind="target"` | Symbol, name, biotype, function, tractability, top associated diseases |
| `genor lookup disease EFO_ID [--n 25]` | `kind="disease"` | Top associated targets with scores |
| `genor lookup variant CHR_POS_REF_ALT [--n 10]` | `kind="variant"` | rsIDs, credible sets containing the variant, top-3 L2G genes per set |
| `genor lookup coloc STUDY_LOCUS_ID [--n 20]` | `kind="coloc"` | Colocalisations of that credible set with other studies (GWAS and QTL) |

The CLI `search` returns 5 hits and has no `--n`; the MCP tool takes `n` (default 10) but no
`entity` filter.

## When to use

- After GWAS loci and candidate genes are known, to see public fine-mapping, L2G, and QTL
  colocalisation for the same variants.
- To check whether a candidate gene is already linked to the disease, and whether it is
  tractable.
- To get disease or target identifiers (EFO/MONDO, Ensembl) for other tools.

## Inputs and prerequisites

- Internet access. No key or account.
- Variant IDs in Open Targets form: `CHR_POS_REF_ALT` on **GRCh38**, e.g. `16_53767042_T_C`
  (rs1421085). See [references/variant-ids.md](references/variant-ids.md) for converting rsIDs and
  GRCh37 positions.
- Ensembl gene IDs (`ENSG...`) for `target`; EFO/MONDO IDs (e.g. `EFO_0004340`, `MONDO_0005148`)
  for `disease`. Get them with `search`.

## Steps

1. Resolve identifiers:
   ```bash
   genor lookup search "body mass index" --entity disease --json
   genor lookup search FTO --entity target --json
   genor lookup search rs1421085 --entity variant --json
   ```
   Search is fuzzy: the first hit for an rsID was the right variant in testing, but other hits were
   unrelated variants. Confirm with `lookup variant` (its `rsIds` field).
2. Variant-level genetics for each lead or credible-set variant:
   ```bash
   genor lookup variant 16_53767042_T_C --n 25 --json
   ```
3. Colocalisations for a credible set of interest (take `studyLocusId` from step 2):
   ```bash
   genor lookup coloc 007688a660da8613b177ce1ad1db36bb --n 50 --json
   ```
4. Target profile for each candidate gene:
   ```bash
   genor lookup target ENSG00000140718 --n 20 --json
   ```
5. Disease-centred view (what Open Targets already ranks for the disease):
   ```bash
   genor lookup disease EFO_0004340 --n 50 --json
   ```
6. Build a triage table (see [references/triage.md](references/triage.md)) and **ask the user**
   which criteria matter for their purpose.
7. Record the query date. Lookups write no manifest, so `genor report build` does not include
   them; add them to the report by hand (see `interpretation-and-reporting`).

## Defaults and why

- `--n` sizes are kept small so results fit an agent's context. Increase them when a variant sits
  in hundreds of credible sets (rs1421085 is in 1,701): the rows returned are the first page, **not
  ranked by p-value or relevance**.
- L2G: GeNor returns the top 3 genes per credible set.
- Tractability: GeNor keeps only tractability entries whose value is true.
- Citation: Buniello A, et al. Open Targets Platform. Nucleic Acids Res 2025;53:D1467-D1475
  (key `opentargets` in GeNor's citation list; add by hand).

## Decision points

- **Ask the user which loci and genes to triage**, and for what purpose (biology vs drug target).
  Triage criteria differ.
- **Ask the user which disease identifier** matches their phenotype when `search` returns several
  (e.g. "major depressive disorder" vs "depressive disorder" vs "unipolar depression").
- **Ask the user before lifting a whole trait to GRCh38** just for lookups; converting a handful of
  variants via `search` is often enough.
- **Ask the user how to weigh tractability vs genetic evidence**; GeNor does not score targets.

## Reading the output

- `variant.credibleSets.count`: how many public credible sets include the variant.
  Each row: `studyLocusId`, `pValueMantissa` x 10^`pValueExponent`, `study.id`,
  `study.traitFromSource`, `study.studyType` (`gwas`, or QTL types such as `eqtl`, `pqtl`, `sqtl`),
  `l2GPredictions.rows` (`score`, `target.approvedSymbol`, `target.id`).
- L2G score: 0-1, a machine-learning prediction that the gene is the causal gene for that credible
  set. A score >= 0.5 is commonly used as a confident assignment. It is specific to a credible
  set, so the same gene can score differently in different studies.
- `credible_set.colocalisation.rows`: `colocalisationMethod`, `h4` (coloc PP.H4), `clpp`
  (eCAVIAR), `otherStudyLocus.study` (`studyType`, `traitFromSource`).
- `target.tractability`: list of `{modality, label}`. Modalities: `SM` small molecule, `AB`
  antibody, `PR` PROTAC, `OC` other clinical modalities. Labels such as "Approved Drug",
  "Advanced Clinical", "High-Quality Pocket", "UniProt loc high conf".
- `associatedDiseases` / `associatedTargets` rows: overall association `score` 0-1 aggregated over
  all evidence types (genetics, somatic, known drugs, pathways, literature, animal models, RNA
  expression). It is a ranking score, not a probability.

Field-by-field details: [references/evidence-types.md](references/evidence-types.md).

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `not_found` (variant) | Variant ID not in Open Targets | Use GRCh38 `CHR_POS_REF_ALT`; try swapping REF/ALT; check with `search`. Variants never in any Open Targets dataset are absent. |
| `not_found` (target / disease) | ID unknown | Use `search` to get the current ID; EFO IDs use an underscore (`EFO_0004340`). |
| `remote_bad_request` | GraphQL error or HTTP 4xx; message included | Usually a malformed ID or an API schema change; report the message. |
| `remote_unavailable` | 5xx or network failure after 3 retries | Retry later; report it to the user. |
| `bad_param` (MCP) | `kind` not one of search, target, disease, variant, coloc | Fix `kind`. |

`lookup coloc` with an unknown `studyLocusId` returns `"credible_set": null` rather than an error.

## What not to conclude

- High L2G does not prove the causal gene; it is a prediction trained on known gene-locus pairs,
  biased towards well-studied genes.
- A high association score does not mean a target is valid; scores rise with literature volume.
- Tractability is about the protein, not about whether modulating it would treat the disease, or
  in which direction.
- Absence from Open Targets is not evidence of no association; coverage depends on which studies
  were ingested and fine-mapped.
- Open Targets credible sets come from public studies that may overlap the user's own GWAS
  samples; agreement is not independent replication.

## Related skills

- `fine-mapping`, `colocalisation`, `gene-pathway-celltype`, `mendelian-randomisation`
  (genetic support for causal direction), `interpretation-and-reporting`, `gwas-concepts`.
- References: [variant-ids.md](references/variant-ids.md),
  [evidence-types.md](references/evidence-types.md), [triage.md](references/triage.md).
