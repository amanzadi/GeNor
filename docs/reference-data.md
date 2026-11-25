# Reference data and external tools

Nothing is downloaded without a command. Methods that need missing data stop with
`refdata_missing` (or `backend_missing`) and the exact command to run; downloads over 100 MB need
`--yes`, so an agent has to ask you first. `genor doctor` shows what each missing piece blocks.

## Reference datasets (`genor refdata list | fetch | remove`)

| Name | Size | Build | Ancestry | Used by | Source / licence |
|---|---|---|---|---|---|
| `ldsc-ldscores-eur` | 13 MB | GRCh37 | EUR | LDSC h2/rg | Zenodo 10515792, CC BY 4.0 |
| `ldsc-weights-eur` | 12 MB | GRCh37 | EUR | LDSC h2/rg | Zenodo 10515792, CC BY 4.0 |
| `g1000-eur` | 288 MB | GRCh37 | EUR | clumping | Zenodo 10515792 (1000 Genomes, open) |
| `g1000-eas` | 323 MB | GRCh37 | EAS | clumping | Zenodo 10515792 (1000 Genomes, open) |
| `chain-hg19-to-hg38` | 0.2 MB | — | — | liftover | UCSC |
| `chain-hg38-to-hg19` | 1.2 MB | — | — | liftover | UCSC |
| `pleiofdr-template` | 262 MB | GRCh37 | EUR | condFDR/conjFDR | precimed |
| `pleiofdr-ref-eur` | 2.4 GB | GRCh37 | EUR | condFDR/conjFDR | precimed |

Every entry has a pinned SHA-256 in `src/genor/refdata/registry.toml`; downloads that do not
match are rejected. The pleioFDR files have no published checksum, so their pins are the hashes
GeNor's maintainers saw on 2026-11-25. The SHA-256 seen is recorded in each run's manifest.

## External tools (`genor backend list | install | add | remove`)

| Name | Mode | Platforms | Used by |
|---|---|---|---|
| `plink2` (alpha 7, 2026-11-07) | downloaded by GeNor, SHA-256 pinned | Linux x86-64 (AVX2 or generic), macOS arm64 / x86-64 | clumping |
| `gcta` (1.95.4) | downloaded, pinned | Linux x86-64, macOS arm64 | (COJO, roadmap) |
| `smr` (1.4.3) | downloaded, pinned | Linux x86-64, macOS arm64 | (SMR, roadmap) |
| `magma` | user installs (licence forbids redistribution) | any | `genor magma` |
| `pleiofdr` | Python package, pinned git commit, Python ≥ 3.13 | any | `genor pleio` |

Resolution order: a path registered with `genor backend add`, then a binary GeNor downloaded,
then `PATH`.
