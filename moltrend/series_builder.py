"""
series_builder.py — Group standardized records into competitor series.

A "series" is a set of compounds from the same company, against the same target,
sharing a common Murcko scaffold. This is the unit a medicinal chemist actually
thinks in — not "building block class" but "Mirati's KRAS G12D series."

Output schema per series (goes into data.json["series"]):
  {
    "id":            "ser_001",
    "target":        "KRAS",
    "company":       "Mirati",
    "company_color": "#00A651",
    "clinical_stage": "Phase 2",          # inferred from source field
    "headline":      "Mirati — KRAS series (fluorinated cyclopropyl)",
    "mechanism":     "Switch-II pocket · gem-difluorocyclopropyl",
    "scaffold_smiles": "FC1(F)CC1",        # Murcko scaffold of majority member
    "compounds": [
      {
        "smiles":    "FC1(CC1)c1ccc(F)cc1",
        "name":      "...",
        "patent_id": "US20230134380",
        "date":      "2023",
        "mw":        228.2,
        "clogp":     2.8,
        "is_clinical": false
      }
    ],
    "patent_ids":   ["US20230134380"],
    "n_compounds":  12,
    "first_seen":   "2022",
    "last_updated": "2024",
    "is_new_this_week": false,
    "tags":         ["fluorinated cyclopropane", "KRAS G12C/D"]
  }
"""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import Optional

logger = logging.getLogger(__name__)

# ── Company color map ─────────────────────────────────────────────────────────
COMPANY_COLORS: dict[str, str] = {
    "Pfizer":                "#003087",
    "Novartis":              "#E74011",
    "AstraZeneca":           "#003865",
    "Roche":                 "#0077C0",
    "Roche / Genentech":     "#0077C0",
    "Merck":                 "#00857C",
    "BMS":                   "#BE2026",
    "Eli Lilly":             "#D52B1E",
    "J&J / Janssen":         "#CA1E24",
    "AbbVie":                "#071D49",
    "Boehringer":            "#DD2323",
    "GSK":                   "#F36633",
    "Sanofi":                "#7A1C8E",
    "Takeda":                "#CC0000",
    "Amgen":                 "#5B4BC8",
    "Gilead":                "#CC3300",
    "Biogen":                "#0066CC",
    "Regeneron":             "#0057A8",
    "Vertex":                "#E6263A",
    "Incyte":                "#00A4E4",
    "Blueprint Medicines":   "#1B9E77",
    "Arvinas":               "#2563EB",
    "Kymera Therapeutics":   "#EA580C",
    "C4 Therapeutics":       "#0891B2",
    "Nurix":                 "#7C3AED",
    "Monte Rosa":            "#D97706",
    "Mirati":                "#1A8C4E",
    "Revolution Medicines":  "#7B2FBE",
    "Araxes Pharma":         "#0F766E",
    "Relay Therapeutics":    "#2563EB",
    "Schrödinger":           "#4F46E5",
}
DEFAULT_COLOR = "#6B7280"


def _company_color(company: str) -> str:
    return COMPANY_COLORS.get(company, DEFAULT_COLOR)


def _infer_clinical_stage(records: list[dict]) -> str:
    """Infer best clinical stage across all records in a series."""
    sources = {r.get("source", "") for r in records}
    if "approval" in sources:
        return "Approved"
    if "clinical" in sources:
        return "Clinical"
    # Check targets for known approved drugs
    targets = set()
    for r in records:
        targets.update(r.get("targets", []))
    return "Preclinical"


def _murcko_scaffold(smiles: str) -> Optional[str]:
    """Return Murcko scaffold SMILES, or None if RDKit unavailable."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold) if scaffold else None
    except Exception:
        return None


def _get_descriptors(smiles: str) -> dict:
    """Return MW and cLogP for a SMILES string."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                "mw":    round(Descriptors.ExactMolWt(mol), 1),
                "clogp": round(Descriptors.MolLogP(mol), 2),
            }
    except Exception:
        pass
    return {"mw": None, "clogp": None}


def _series_headline(company: str, target: str, tags: list[str]) -> str:
    parts = [f"{company} — {target} series"]
    if tags:
        parts.append(f"({tags[0]})")
    return " ".join(parts)


def _series_mechanism(tags: list[str], records: list[dict]) -> str:
    """Build a short mechanism string from tags and class info."""
    # Infer from class names if present
    classes = set()
    for r in records:
        for c in r.get("classes", []):
            classes.add(c)

    parts = []
    if "Acrylamide warhead" in classes or "Cyanoacrylamide warhead" in classes:
        parts.append("Covalent irreversible")
    if "Sulfonyl fluoride (SuFEx)" in classes:
        parts.append("SuFEx covalent")
    if "Cereblon glutarimide" in classes:
        parts.append("CRBN binder")
    if "PROTAC linker" in classes or "Heterobifunctional linker" in classes:
        parts.append("Bifunctional degrader")
    if "BCP bioisostere" in classes:
        parts.append("BCP sp3 bioisostere")
    if "Bicyclo[1.1.0]butane (BCB)" in classes:
        parts.append("BCB strain-release warhead")
    if tags:
        parts.append(tags[0])
    return " · ".join(parts[:2]) if parts else "Small molecule"


def build_series(records: list[dict], min_per_series: int = 2) -> list[dict]:
    """
    Group records into competitor series.

    Groups by: (target, assignee) → Murcko scaffold cluster
    Minimum min_per_series compounds required to surface a series.

    Returns list of series dicts for data.json.
    """
    # Filter to records that have both target and assignee
    attributed = [
        r for r in records
        if r.get("assignee") and r.get("targets")
    ]

    if not attributed:
        logger.debug("series_builder: no attributed records to group")
        return []

    # ── Phase 1: group by (primary_target, assignee) ─────────────────────────
    buckets: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in attributed:
        company = r["assignee"]
        # Use first target as primary key
        target  = r["targets"][0] if r["targets"] else "Unknown"
        buckets[(target, company)].append(r)

    # ── Phase 2: within each bucket, optionally sub-cluster by Murcko scaffold
    series_list: list[dict] = []
    ser_idx = 1

    for (target, company), recs in sorted(buckets.items()):
        if len(recs) < min_per_series:
            continue

        # Compute scaffold for each record
        scaffold_buckets: dict[str, list[dict]] = defaultdict(list)
        for r in recs:
            scaffold = _murcko_scaffold(r.get("smiles", "")) or "__none__"
            scaffold_buckets[scaffold].append(r)

        # Keep only scaffold sub-groups with enough members
        sub_groups = [
            (scaffold, sub_recs)
            for scaffold, sub_recs in scaffold_buckets.items()
            if len(sub_recs) >= min_per_series or scaffold == "__none__"
        ]

        # If all records have no scaffold (e.g. fragments), treat as one group
        if not sub_groups:
            sub_groups = [("__none__", recs)]

        for scaffold_smi, sub_recs in sub_groups:
            if len(sub_recs) < min_per_series:
                continue

            # Build compound list
            compounds = []
            seen_smi = set()
            for r in sorted(sub_recs, key=lambda x: x.get("date", "0"), reverse=True):
                smi = r.get("smiles", "")
                if smi in seen_smi:
                    continue
                seen_smi.add(smi)
                desc = _get_descriptors(smi)
                compounds.append({
                    "smiles":      smi,
                    "name":        r.get("name") or None,
                    "patent_id":   r.get("ref_id") or None,
                    "date":        r.get("date") or None,
                    "mw":          desc["mw"],
                    "clogp":       desc["clogp"],
                    "is_clinical": r.get("source") in ("approval", "clinical"),
                })

            # Tags from classes
            tags = list({
                c for r in sub_recs
                for c in r.get("classes", [])
                if c not in ("Spirocyclic scaffold",)
            })[:3]

            # Collect patent IDs
            patent_ids = sorted({
                r.get("ref_id", "") for r in sub_recs if r.get("ref_id")
            })

            dates = sorted([r.get("date", "") for r in sub_recs if r.get("date")])

            series_list.append({
                "id":             f"ser_{ser_idx:03d}",
                "target":         target,
                "company":        company,
                "company_color":  _company_color(company),
                "clinical_stage": _infer_clinical_stage(sub_recs),
                "headline":       _series_headline(company, target, tags),
                "mechanism":      _series_mechanism(tags, sub_recs),
                "scaffold_smiles": scaffold_smi if scaffold_smi != "__none__" else (compounds[0]["smiles"] if compounds else ""),
                "compounds":      compounds[:30],  # cap for JSON size
                "patent_ids":     patent_ids[:5],
                "n_compounds":    len(compounds),
                "first_seen":     dates[0][:4] if dates else None,
                "last_updated":   dates[-1][:4] if dates else None,
                "is_new_this_week": False,  # pipeline sets this via date comparison
                "tags":           tags,
            })
            ser_idx += 1

    # Sort: clinical first, then by company name
    stage_order = {"Approved": 0, "Clinical": 1, "Preclinical": 2}
    series_list.sort(key=lambda s: (stage_order.get(s["clinical_stage"], 3), s["company"]))

    logger.info("series_builder: %d series from %d attributed records", len(series_list), len(attributed))
    return series_list
