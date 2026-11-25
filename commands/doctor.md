---
description: Check GeNor's setup (backends, reference data) and say what is blocked
allowed-tools: Bash(genor *)
---
Run `genor doctor --json` and summarise for the user in a few lines: platform, which analyses work
now, and each missing backend or reference dataset with what it blocks and its size. Do not
download or install anything; ask first.
