# Overleaf figure sync

`tools/sync_overleaf_figures.py` copies only configured figure files from the
monorepo into the proceeding's `figures/` directory. It compares SHA-256
hashes, preserves unchanged files, never deletes destination files, and writes
`figures/SOURCES.md` with each source, result run and date (when applicable),
hash, repository commit, and last content-changing sync time.

The mapping lives in `tools/overleaf_figures.json`:

```json
{
  "figdir": "figures",
  "mappings": [
    {
      "dst": "figures/example.png",
      "src": "results/inversion/2026-05-09_example/example.png",
      "run": "results/inversion/2026-05-09_example"
    }
  ],
  "unresolved": []
}
```

Paths beginning with `results/` are resolved through `data/`; other source
paths are relative to the repository root. Use `run: null` for tracked files
that are not result-run outputs. The `unresolved` list records figures with
zero or multiple monorepo candidates for manual review.

The destination defaults to the Dropbox Overleaf project. Set
`MUOGRAPHY_OVERLEAF_DIR` to use another project directory.

- `make sync-figures-dry` previews copies and the provenance file.
- `make sync-figures` performs the configured one-way sync.
- `make check-figures` checks references and whether mapped files are current.

Operating rule: close Overleaf, run `make sync-figures`, wait for Dropbox to
finish syncing, then reopen Overleaf. Never run the sync while editing in
Overleaf. Review unresolved figures before relying on `make check-figures`.
