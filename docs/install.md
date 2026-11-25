# Installation

GeNor needs only [uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

uv fetches a suitable Python and GeNor's dependencies on first use. Supported platforms: Linux and
macOS (x86-64 and Apple silicon). On Windows, install uv inside WSL and follow the Linux
instructions.

Until GeNor is on PyPI, replace `genor` in `uvx genor ...` with
`--from git+https://github.com/amanzadi/Genor genor`, e.g.
`uvx --from git+https://github.com/amanzadi/Genor genor doctor`.

Claude Code and Claude Desktop are the fully tested surfaces. Codex and opencode use the same MCP
server and skills; they are smoke-tested only.

## Claude Code (plugin)

```text
/plugin marketplace add amanzadi/Genor
/plugin install genor@genor
```

or from a shell: `claude plugin marketplace add amanzadi/Genor && claude plugin install genor@genor`.

The plugin provides:

- the 12 GeNor skills (auto-loaded when relevant; start with `genor`)
- the `genor` MCP server (`uvx --from ${CLAUDE_PLUGIN_ROOT} genor mcp`)
- `genor` on the PATH of Claude's Bash tool (`bin/genor`, a thin `uvx` wrapper)
- slash commands: `/genor:doctor`, `/genor:pipeline`, `/genor:mr`, `/genor:rg`, `/genor:example`

The first call builds GeNor's environment (about a minute); later calls start in seconds.

Without the plugin, use the CLI and copy the skills into the project:

```bash
uv tool install git+https://github.com/amanzadi/Genor     # puts `genor` on PATH
genor skills install --agent claude                       # -> .claude/skills (or --scope user)
claude mcp add genor -- genor mcp                         # optional: MCP tools too
```

## Claude Desktop

Claude Desktop reaches GeNor through MCP only. Edit `claude_desktop_config.json`
(Settings → Developer → Edit Config):

```json
{
  "mcpServers": {
    "genor": {
      "command": "uvx",
      "args": ["--from", "git+https://github.com/amanzadi/Genor", "genor", "mcp"],
      "env": { "GENOR_WORKDIR": "/Users/you/genor-workspace" }
    }
  }
}
```

Claude Desktop does not inherit your shell PATH: if `uvx` is not found, use its absolute path
(`which uvx`, typically `~/.local/bin/uvx`). `GENOR_WORKDIR` is where runs are written and where
relative file paths are resolved; put your summary statistics there or give absolute paths.

Desktop chats have no skills, so the server tells the model to call `get_guidance("genor")`
first; the skill text is served over MCP.

An experimental Desktop extension manifest is in `packaging/mcpb/` (see its README).

## Codex (CLI and desktop app)

`~/.codex/config.toml`:

```toml
[mcp_servers.genor]
command = "uvx"
args = ["--from", "git+https://github.com/amanzadi/Genor", "genor", "mcp"]
startup_timeout_sec = 120   # first start builds the environment
tool_timeout_sec = 60       # GeNor tools return within ~50 s and hand back job ids

[mcp_servers.genor.env]
GENOR_WORKDIR = "/home/you/genor-workspace"
```

Skills (Codex reads `.agents/skills` in the repo and `~/.agents/skills`):

```bash
uvx --from git+https://github.com/amanzadi/Genor genor skills install --agent codex            # project
uvx --from git+https://github.com/amanzadi/Genor genor skills install --agent codex --scope user
```

If the Codex desktop app does not load skills, the MCP `get_guidance` tool serves the same text.

## opencode

`opencode.json` (project) or `~/.config/opencode/opencode.json`. Recent opencode versions:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "servers": {
      "genor": {
        "type": "local",
        "command": ["uvx", "--from", "git+https://github.com/amanzadi/Genor", "genor", "mcp"]
      }
    }
  }
}
```

Older opencode versions put the server directly under `"mcp"` with `"enabled": true`; check
`opencode --version` and its docs. Skills: `genor skills install --agent opencode` writes to
`.agents/skills`, which opencode also reads.

## Plain CLI

```bash
uvx genor doctor                 # one-off
uv tool install genor            # persistent `genor` command
genor --help
```

## Optional pieces

| Need | Command | Notes |
|---|---|---|
| Clumping / LD | `genor backend install plink2` and `genor refdata fetch g1000-eur --yes` | 288 MB panel |
| LDSC | `genor refdata fetch ldsc-ldscores-eur ldsc-weights-eur` | 25 MB |
| condFDR / conjFDR | run GeNor on Python 3.13 with the pleioFDR port (below) | 262 MB template + 2.4 GB reference, ~12 GB RAM |
| MAGMA | download from https://cncr.nl/research/magma/, then `genor backend add magma /path/to/magma` | licence forbids bundling |
| COJO, SMR binaries | `genor backend install gcta` / `smr` | commands that drive them are on the roadmap |

pleioFDR (Python port, pinned commit):

```bash
uv tool install --python 3.13 git+https://github.com/amanzadi/Genor \
  --with "pleiofdr @ git+https://github.com/espenhgn/pleiofdr@45e03b9c98c48ce9ef1496c4173a5581b23a853d"
```

## Where things are stored

| What | Default | Override |
|---|---|---|
| Reference data, downloaded binaries | `~/Library/Application Support/genor` (macOS), `~/.local/share/genor` (Linux) | `GENOR_DATA_DIR` or `GENOR_HOME` |
| Backend registrations | platform config dir | `GENOR_CONFIG_DIR` |
| Background jobs | platform state dir | `GENOR_JOBS_DIR` |
| Run directories | `./genor_runs/` or `--out` | — |

On a cluster, set `GENOR_HOME=/shared/project/genor` so reference data is downloaded once.
