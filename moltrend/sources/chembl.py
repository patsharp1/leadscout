"""
chembl.py — ChEMBL source adapter.

Pulls two layers:
  1. Approval/clinical layer: molecules with max_phase >= 1.
     source tag: "approval" (phase 4) / "clinical" (phase 1-3)
  2. Literature/bioactivity layer: recent bioactivity records with SMILES +
     year from linked publication.
     source tag: "paper"

No API key required. Uses the ChEMBL REST API:
  https://www.ebi.ac.uk/chembl/api/data/

Caching: parquet files in .moltrend_cache/ to make reruns cheap.
"""
from __future__ import annotations

import json
import logging
import os
import time
import urllib.request
import urllib.parse
from datetime import datetime

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_DIR   = ".moltrend_cache"
PAGE_SIZE   = 1000    # records per API page

# Literature: restrict to recent years and drug-discovery relevant assay types
LITERATURE_MIN_YEAR = 2017
ACTIVITY_TYPES      = ["IC50", "Ki", "Kd", "EC50", "GI50"]
MAX_ACTIVITY_PAGES  = 20   # ~20k activity records — manageable, cache aggressively


def fetch() -> list[dict]:
    """Fetch both layers from ChEMBL. Returns list of record dicts."""
    records = []
    records.extend(_fetch_molecules())
    records.extend(_fetch_activities())
    return records


# ── Layer 1: Molecules (approval / clinical) ──────────────────────────────────

def _fetch_molecules() -> list[dict]:
    cache_path = _cache_path("chembl_molecules.json")
    if os.path.exists(cache_path):
        logger.info("  [cache] %s", cache_path)
        with open(cache_path) as f:
            return json.load(f)

    logger.info("  Fetching ChEMBL molecules (max_phase>=1)…")
    records = []
    offset  = 0

    while True:
        url = (f"{CHEMBL_BASE}/molecule.json"
               f"?max_phase__gte=1"
               f"&limit={PAGE_SIZE}&offset={offset}"
               f"&format=json")
        data = _get_json(url)
        if data is None:
            break
        mols = data.get("molecules", [])
        for m in mols:
            r = _molecule_to_record(m)
            if r:
                records.append(r)

        total = data.get("page_meta", {}).get("total_count", 0)
        offset += PAGE_SIZE
        if offset >= total:
            break
        time.sleep(0.2)  # polite rate-limiting

    logger.info("  ChEMBL molecules: %d records", len(records))
    _write_cache(cache_path, records)
    return records


def _molecule_to_record(m: dict) -> dict | None:
    smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
    if not smi:
        return None
    phase = m.get("max_phase", 0)
    year  = m.get("first_approval")
    return {
        "smiles":   smi,
        "source":   "approval" if phase == 4 else "clinical",
        "date":     str(year) if year else str(datetime.now().year - 5),
        "name":     m.get("pref_name") or m.get("molecule_chembl_id"),
        "ref_id":   m.get("molecule_chembl_id"),
        "assignee": None,
        "targets":  [],  # enriched lazily — expensive per molecule
    }


# ── Layer 2: Bioactivity records (literature) ─────────────────────────────────

def _fetch_activities() -> list[dict]:
    cache_path = _cache_path("chembl_activities.json")
    if os.path.exists(cache_path):
        logger.info("  [cache] %s", cache_path)
        with open(cache_path) as f:
            return json.load(f)

    logger.info("  Fetching ChEMBL activities (year>=%d)…", LITERATURE_MIN_YEAR)
    records = []
    offset  = 0
    pages   = 0

    activity_filter = "&".join([f"standard_type={t}" for t in ACTIVITY_TYPES])

    while pages < MAX_ACTIVITY_PAGES:
        url = (f"{CHEMBL_BASE}/activity.json"
               f"?{activity_filter}"
               f"&document_year__gte={LITERATURE_MIN_YEAR}"
               f"&limit={PAGE_SIZE}&offset={offset}&format=json")
        data = _get_json(url)
        if data is None:
            break
        acts = data.get("activities", [])
        for a in acts:
            r = _activity_to_record(a)
            if r:
                records.append(r)
        total = data.get("page_meta", {}).get("total_count", 0)
        offset += PAGE_SIZE
        pages  += 1
        if offset >= total:
            break
        time.sleep(0.3)

    logger.info("  ChEMBL activities: %d records", len(records))
    _write_cache(cache_path, records)
    return records


def _activity_to_record(a: dict) -> dict | None:
    smi = a.get("canonical_smiles")
    if not smi:
        return None
    year = a.get("document_year")
    target = a.get("target_pref_name") or ""
    return {
        "smiles":   smi,
        "source":   "paper",
        "date":     str(year) if year else str(LITERATURE_MIN_YEAR),
        "name":     a.get("molecule_pref_name"),
        "ref_id":   a.get("molecule_chembl_id"),
        "assignee": None,
        "targets":  [target] if target else [],
    }


# ── Utilities ─────────────────────────────────────────────────────────────────

def _get_json(url: str, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            logger.debug("GET failed (attempt %d): %s — %s", attempt + 1, url[:80], e)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def _cache_path(name: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def _write_cache(path: str, data) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
