# Open Targets variant IDs

Open Targets identifies variants as `CHR_POS_REF_ALT`:

- `CHR`: 1-22, X, Y (no `chr` prefix).
- `POS`: 1-based position on **GRCh38**.
- `REF`: the GRCh38 reference allele; `ALT`: the alternative allele.

Example: rs1421085 is `16_53767042_T_C` (GRCh37 position 16:53800954).

GWAS files use `effect_allele` / `other_allele`, which are not REF/ALT. The effect allele may be
REF or ALT. When building an ID from a GWAS row, try `CHR_POS_OA_EA` and `CHR_POS_EA_OA`; only one
exists.

## Getting from what you have to an Open Targets ID

### From an rsID (fastest; no liftover)

```bash
genor lookup search rs1421085 --entity variant --json
```

Take the hit whose ID then returns the same rsID in `genor lookup variant <id>` (`rsIds` field).
Do not take the first hit blindly: search is fuzzy and returned unrelated variants after the
correct one in testing.

Alternatively, `genor lookup vep rs1421085 --build GRCh38` gives the GRCh38 `location`
(chromosome:position) but not REF/ALT.

### From GRCh37 positions (whole trait)

Lift the trait to GRCh38 (needs the chain file, about 0.2 MB):

```bash
genor refdata fetch chain-hg19-to-hg38
genor sumstats liftover runs/bmi --to GRCh38 --out runs/bmi_GRCh38 --json
```

Check `n_unmapped`; the warning `liftover_many_unmapped` appears if more than 5% fail. Liftover
changes positions only, not alleles; strand flips between builds are rare but exist. Ask the user
before lifting a large file just for a few lookups.

## Multi-allelic sites and indels

- Multi-allelic sites have several IDs at one position (e.g. `1_100_A_G` and `1_100_A_T`).
- Indels use the VCF convention with an anchor base (e.g. `1_100_A_AT`). GWAS files sometimes
  encode indels as `I`/`D` or `R`/`D`; these cannot be mapped without the reference sequence.

## When a variant is `not_found`

- Wrong build (GRCh37 position used).
- REF/ALT reversed.
- Variant not in any Open Targets dataset (rare, or never in a fine-mapped credible set or the
  variant index). This is not evidence of no association.
