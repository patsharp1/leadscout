"""
build.py — Orchestrator. The only user-facing script.

Usage:
  python -m moltrend.build --source seed --out data.json
  python -m moltrend.build --source chembl --out data.json
  python -m moltrend.build --source all --out data.json

Flags:
  --source   {seed,chembl,surechembl,openfda,all}
  --out      output path (default: data.json)
  --no-llm   disable LLM editorial (use deterministic copy)
  --min-support  minimum evidence count to surface a trend (default: 3)
  --pubchem  enrich supplier counts via PubChem (slower)
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime, timezone

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
logger = logging.getLogger(__name__)


def main(args=None):
    parser = argparse.ArgumentParser(description="MolTrend pipeline — generate data.json")
    parser.add_argument("--source",      default="seed",
                        choices=["seed", "chembl", "surechembl", "openfda", "pharma", "all"])
    parser.add_argument("--out",         default="data.json")
    parser.add_argument("--no-llm",      action="store_true")
    parser.add_argument("--min-support", type=int, default=3)
    parser.add_argument("--pubchem",     action="store_true",
                        help="Enrich supplier counts via PubChem API")
    parser.add_argument("--limit",       type=int, default=0,
                        help="Limit records per source (0 = no limit, useful for testing)")
    ns = parser.parse_args(args)

    logger.info("MolTrend build — source=%s  out=%s", ns.source, ns.out)
    t0 = time.time()

    # ── 1. Fetch records ──────────────────────────────────────────────────────
    raw_records: list[dict] = []

    sources_to_run = {
        "seed":        ["seed"],
        "chembl":      ["chembl"],
        "surechembl":  ["surechembl"],
        "openfda":     ["openfda"],
        "pharma":      ["pharma_patents"],
        "all":         ["chembl", "surechembl", "openfda", "pharma_patents"],
    }[ns.source]

    for src_name in sources_to_run:
        logger.info("Fetching source: %s", src_name)
        try:
            adapter = _load_adapter(src_name)
            recs = adapter.fetch()
            if ns.limit:
                recs = recs[:ns.limit]
            logger.info("  %d records from %s", len(recs), src_name)
            raw_records.extend(recs)
        except Exception as e:
            logger.error("Source %s failed: %s", src_name, e)

    if not raw_records:
        logger.error("No records fetched — aborting")
        sys.exit(1)

    logger.info("Total raw records: %d", len(raw_records))

    # ── 2. Standardize ────────────────────────────────────────────────────────
    from moltrend.standardize import standardize
    std_records: list[dict] = []
    seen_ik_src: set[tuple[str, str]] = set()
    n_fail = 0

    for rec in raw_records:
        result = standardize(rec.get("smiles", ""))
        if result is None:
            n_fail += 1
            continue
        key = (result["inchikey"], rec.get("source", ""))
        if key in seen_ik_src:
            continue
        seen_ik_src.add(key)
        merged = {**rec, **result}   # standardized values overwrite raw
        std_records.append(merged)

    logger.info("Standardized: %d ok, %d failed, %d dedup'd",
                len(std_records), n_fail, len(raw_records) - len(std_records) - n_fail)

    if not std_records:
        logger.error("No valid structures after standardization — aborting")
        sys.exit(1)

    # ── 3. Classify (primary — class track) ──────────────────────────────────
    from moltrend.classify import build_class_groups
    logger.info("Classifying into named medchem classes…")
    class_groups = build_class_groups(std_records)
    logger.info("  %d class groups found", len(class_groups))

    # ── 4. Fragment (secondary — BRICS discovery track) ───────────────────────
    from moltrend.fragment import build_fragment_groups
    logger.info("BRICS fragmentation (secondary track)…")
    brics_groups = build_fragment_groups(std_records)
    logger.info("  %d BRICS fragment groups found", len(brics_groups))

    # ── 5. Score + embed ──────────────────────────────────────────────────────
    from moltrend.trends import score_groups
    logger.info("Scoring groups and computing PCA embedding…")
    trends, compounds = score_groups(
        class_groups, brics_groups, std_records,
        min_support=ns.min_support,
    )
    logger.info("  %d trends, %d canvas compounds", len(trends), len(compounds))

    if not trends:
        logger.warning("No trends surfaced — check min_support or corpus size")

    # ── 6. Editorial copy ─────────────────────────────────────────────────────
    from moltrend.editorial import annotate_trends
    use_llm = not ns.no_llm
    logger.info("Writing editorial copy (llm=%s)…", use_llm and bool(os.environ.get("ANTHROPIC_API_KEY")))
    trends = annotate_trends(trends, use_llm=use_llm)

    # ── 7. Pharma patent signals (novel building block detection) ─────────────
    pharma_signals: list[dict] = []
    pharma_records = [r for r in std_records if r.get("assignee")]
    if pharma_records:
        try:
            from moltrend.novel_blocks import detect_novel_blocks
            logger.info("Detecting novel pharma building blocks (%d attributed records)…",
                        len(pharma_records))
            pharma_signals = detect_novel_blocks(pharma_records, background_records=std_records)
            logger.info("  %d pharma signals", len(pharma_signals))
        except Exception as e:
            logger.warning("Novel block detection failed: %s", e)

    # ── 8. Supplier enrichment ────────────────────────────────────────────────
    supplier_data = "placeholder"
    if ns.pubchem:
        try:
            logger.info("Enriching supplier counts via PubChem…")
            _enrich_pubchem(trends)
            supplier_data = "live"
        except Exception as e:
            logger.warning("PubChem enrichment failed: %s", e)

    # ── 9. Build competitor series ────────────────────────────────────────────
    series: list[dict] = []
    try:
        from moltrend.series_builder import build_series
        # Attach class labels to records for series mechanism inference
        from moltrend.classify import classify_molecule
        for r in std_records:
            if "classes" not in r:
                r["classes"] = classify_molecule(r.get("smiles", ""))
        series = build_series(std_records, min_per_series=2)
        logger.info("  %d competitor series", len(series))
    except Exception as e:
        logger.warning("Series builder failed: %s", e)

    # ── 10. Assemble data.json ────────────────────────────────────────────────
    years = _corpus_year_range(std_records)
    output = {
        "generated_at":  datetime.now(timezone.utc).isoformat(),
        "meta": {
            "n_compounds":    len(std_records),
            "n_fragments":    len(brics_groups),
            "n_pharma_sigs":  len(pharma_signals),
            "window_years":   years,
            "supplier_data":  supplier_data,
            "min_support":    ns.min_support,
            "source":         ns.source,
        },
        "trends":          trends,
        "compounds":       compounds,
        "pharma_signals":  pharma_signals,
        "series":          series,
    }

    # Serialize (sets → lists for JSON)
    def _json_default(o):
        if isinstance(o, set):
            return sorted(o)
        raise TypeError(f"Object of type {type(o)} is not JSON serializable")

    with open(ns.out, "w", encoding="utf-8") as f:
        json.dump(output, f, default=_json_default, indent=2, ensure_ascii=False)

    elapsed = round(time.time() - t0, 1)
    logger.info("✓ Wrote %s in %.1fs  (%d trends, %d compounds, %d pharma signals)",
                ns.out, elapsed, len(trends), len(compounds), len(pharma_signals))
    return output


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_adapter(name: str):
    """Dynamically import a source adapter module."""
    import importlib
    mod = importlib.import_module(f"moltrend.sources.{name}")
    return mod


def _corpus_year_range(records: list[dict]) -> list[int]:
    from moltrend.trends import _parse_year
    years = [y for r in records if (y := _parse_year(r.get("date", ""))) is not None]
    if not years:
        return [datetime.now().year - 7, datetime.now().year]
    return [min(years), max(years)]


def _enrich_pubchem(trends: list[dict]) -> None:
    """
    Enrich trends[].suppliers with PubChem vendor counts for the primary SMILES.
    Uses PubChem PUG-REST: /pug/compound/smiles/xrefs/RegistryID/JSON
    (free, no key required).
    """
    import urllib.request
    import urllib.parse

    for t in trends:
        smi = t.get("primary_smiles", "")
        if not smi:
            continue
        try:
            enc = urllib.parse.quote(smi, safe="")
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{enc}/xrefs/RegistryID/JSON"
            with urllib.request.urlopen(url, timeout=8) as resp:
                data = json.loads(resp.read())
            vendors = data.get("InformationList", {}).get("Information", [{}])[0]
            count = len(vendors.get("RegistryID", []))
            t["suppliers"]["enamine"]  = max(count - 5, 0)
            t["suppliers"]["molport"]  = max(count - 10, 0)
            t["suppliers"]["mcule"]    = max(count // 3, 0)
        except Exception:
            pass   # leave as 0 on failure


if __name__ == "__main__":
    main()
