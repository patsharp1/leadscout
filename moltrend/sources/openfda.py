"""
openfda.py — openFDA freshness overlay.

Pulls the most recently approved NMEs from the FDA Drugs@FDA dataset.
Structures are thin in openFDA — we cross-reference compound names to
ChEMBL/PubChem for canonical SMILES.

Role: freshness check for the last 12 months, not a primary trend body.
Source tag: "approval"

API: https://api.fda.gov/drug/drugsfda.json
No key required (rate limit: 240/min unauthenticated).
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

OPENFDA_URL  = "https://api.fda.gov/drug/drugsfda.json"
PUBCHEM_URL  = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
CACHE_DIR    = ".moltrend_cache"
LOOKBACK_YEARS = 2    # fetch approvals from last N years
MAX_APPROVALS  = 100  # NMEs per year is ~50; 2 years = ~100


def fetch() -> list[dict]:
    """Fetch recent NME approvals with SMILES from PubChem cross-reference."""
    cache_path = _cache_path("openfda_records.json")
    if os.path.exists(cache_path):
        logger.info("  [cache] %s", cache_path)
        with open(cache_path) as f:
            return json.load(f)

    nmes  = _fetch_nmes()
    recs  = []
    names_seen: set[str] = set()

    for nme in nmes:
        name = nme.get("name", "")
        if not name or name.lower() in names_seen:
            continue
        names_seen.add(name.lower())
        smi = _name_to_smiles_pubchem(name)
        if not smi:
            continue
        recs.append({
            "smiles":   smi,
            "source":   "approval",
            "date":     nme.get("date", str(datetime.now().year)),
            "name":     name,
            "ref_id":   nme.get("application_number", ""),
            "assignee": nme.get("sponsor", None),
            "targets":  [],
        })
        time.sleep(0.1)   # polite rate-limiting toward PubChem

    logger.info("openFDA: %d NMEs with SMILES", len(recs))
    _write_cache(cache_path, recs)
    return recs


def _fetch_nmes() -> list[dict]:
    """Fetch recent NDA/BLA approvals from openFDA."""
    year_cutoff = datetime.now().year - LOOKBACK_YEARS
    date_str    = f"{year_cutoff}0101"

    url = (f"{OPENFDA_URL}"
           f"?search=submissions.submission_type:NDA+AND+submissions.submission_status_date:[{date_str}+TO+99991231]"
           f"&limit={MAX_APPROVALS}")
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        logger.warning("openFDA fetch failed: %s", e)
        return []

    nmes = []
    for result in data.get("results", []):
        # Extract active ingredient name
        products = result.get("products", [])
        name = ""
        for p in products:
            ing = p.get("active_ingredients", [])
            if ing:
                name = ing[0].get("name", "")
                break
        if not name:
            continue

        # Extract approval date
        sponsor = result.get("sponsor_name", "")
        submissions = result.get("submissions", [])
        date = str(datetime.now().year)
        for sub in submissions:
            if sub.get("submission_status") == "AP":
                d = sub.get("submission_status_date", "")
                if d:
                    date = d[:4]
                    break

        nmes.append({
            "name":    name,
            "date":    date,
            "sponsor": sponsor,
            "application_number": result.get("application_number", ""),
        })

    return nmes


def _name_to_smiles_pubchem(drug_name: str) -> str | None:
    """Cross-reference a drug name to canonical SMILES via PubChem."""
    enc = urllib.parse.quote(drug_name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc}/property/CanonicalSMILES/JSON"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("CanonicalSMILES")
    except Exception:
        pass
    return None


def _cache_path(name: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def _write_cache(path: str, data) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
