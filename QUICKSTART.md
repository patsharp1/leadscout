# MolTrend — Quickstart

## Local dev

```bash
pip install rdkit anthropic

# Fastest: seed corpus (offline, no network)
python -m moltrend.build --source seed --out data.json

# ChEMBL only (fetches ~5k approved + clinical structures)
python -m moltrend.build --source chembl --out data.json

# Full pipeline (ChEMBL + SureChEMBL patents + openFDA overlay)
python -m moltrend.build --source all --out data.json

# With LLM editorial copy and PubChem supplier counts
ANTHROPIC_API_KEY=sk-... python -m moltrend.build --source all --out data.json --pubchem
```

Then open `moltrend.html` in a browser. It fetches `data.json` from the same directory.

> **Note:** Browsers block `fetch()` for local files. Use a simple server:
> ```bash
> python -m http.server 8080
> # open http://localhost:8080/moltrend.html
> ```

## CLI flags

| Flag | Default | Description |
|------|---------|-------------|
| `--source` | `seed` | `seed`, `chembl`, `surechembl`, `openfda`, `all` |
| `--out` | `data.json` | Output path |
| `--min-support` | `3` | Minimum molecules to surface a trend |
| `--no-llm` | off | Disable LLM editorial (always use templates) |
| `--pubchem` | off | Enrich supplier counts via PubChem API |
| `--limit N` | 0 | Limit records per source (0 = no limit; useful for testing) |

## Deploy (GitHub Pages / Vercel / Netlify)

1. Push this repo to GitHub.
2. Enable GitHub Pages from the `main` branch root.
3. Add secrets: `ANTHROPIC_API_KEY` (optional), `MCULE_API_KEY` (optional).
4. The GitHub Actions workflow at `.github/workflows/update.yml` runs nightly at 06:00 UTC, regenerates `data.json`, and commits it. GitHub Pages auto-redeploys.

## Extending the taxonomy

Edit `moltrend/classify.py` → `CLASS_SMARTS` dict. Each entry is a display name → list of SMARTS patterns. The class appears in the feed as soon as `min_support` molecules match it.

```python
CLASS_SMARTS["Macrocyclic peptide"] = ["C1NC(=O)CCCCCCC1"]
```

## Adding a new source

Create `moltrend/sources/mydb.py` with a single `fetch() -> list[dict]` function. Each record must follow the schema:

```python
{
    "smiles":   str,                          # raw; standardized downstream
    "source":   "patent|paper|approval|clinical",
    "date":     "YYYY-MM-DD" | "YYYY",
    "assignee": str | None,
    "targets":  [str],
    "ref_id":   str,
    "name":     str | None,
}
```

Then run: `python -m moltrend.build --source mydb`
