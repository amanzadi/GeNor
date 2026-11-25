# Claude Desktop extension (experimental)

`manifest.json` describes a `.mcpb` bundle that launches `uvx ... genor mcp`. The MCPB `uv`
server type is new and this manifest has not been tested in Claude Desktop yet; the
JSON config in `docs/install.md` is the supported route.

Build: `npx @anthropic-ai/mcpb pack packaging/mcpb genor.mcpb`, then open `genor.mcpb` in Claude
Desktop. Requires `uv` on the machine.
