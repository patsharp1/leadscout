"""
novel_blocks.py — Detect novel building blocks in pharma patent records.

For each pharma-sourced record, BRICS-decomposes the molecule and
tags each fragment with:
  - which company filed it
  - which targets the patent covers
  - how many patents it appears in (support)
  - a novelty score (inverse frequency in non-pharma corpus)

Outputs a `pharma_signals` list for data.json:

[
  {
    "id": "ps1",
    "fragment_smiles": "FC1CC1",
    "bb_class": "Fluorinated cyclopropane",
    "assignees": ["J&J / Janssen", "Roche"],
    "targets":   ["KRAS"],
    "support":   4,
    "headline":  "Fluorinated cyclopropane in J&J filing targeting KRAS",
    "badge":     "novel",
    "first_seen": "2023",
    "sparkline": [0,0,1,3,4],
  },
  ...
]

Integration: called from build.py after the normal pipeline, merged into data.json.
"""
from __future__ import annotations

import logging
import re
from collections import defaultdict
from typing import Optional

logger = logging.getLogger(__name__)

# ── SMARTS for named building-block classes ───────────────────────────────────
# These augment classify.py with fragment-level patterns for the pharma feed.
# Each entry: (display_name, list_of_SMARTS, priority)
# Priority: higher = shown first if multiple classes match

BB_SMARTS: list[tuple[str, list[str], int]] = [
    # Strained / fluorinated carbocycles (the "rising" story right now)
    ("Fluorinated cyclopropane",  ["FC1CC1", "C(F)(F)C1CC1", "[F]C1CC1"],                  10),
    ("gem-Difluorocyclopropane",  ["FC1(F)CC1"],                                             12),
    ("Fluorinated cyclobutane",   ["FC1CCC1", "C(F)(F)C1CCC1"],                              9),
    ("gem-Difluorocyclobutane",   ["FC1(F)CCC1"],                                            11),
    ("Spiro[2.2]pentane",         ["C1CC12CC2"],                                              9),
    ("Bicyclo[1.1.0]butane (BCB)", ["C12CC1C2"],                                             10),
    ("Bicyclo[1.1.1]pentane (BCP)", ["C12CCC1C2"],                                           11),
    ("Oxetane",                   ["C1COC1"],                                                  7),
    ("Azetidine",                 ["C1CNC1"],                                                  7),
    # Warheads
    ("Acrylamide warhead",        ["[CX3]=[CX3]C(=O)[NX3]", "C#CC(=O)[NX3]"],              10),
    ("Cyanoacrylamide warhead",   ["N#CC=CC(=O)[NX3]"],                                     11),
    ("Sulfonyl fluoride (SuFEx)", ["S(=O)(=O)F"],                                            10),
    ("Chloroacetamide",           ["ClCC(=O)[NX3]"],                                          9),
    ("Epoxide warhead",           ["C1OC1"],                                                   8),
    # E3 ligase handles
    ("Glutarimide (CRBN)",        ["O=C1CCCC(=O)N1"],                                       10),
    ("Hydroxamic acid",           ["C(=O)NO"],                                                8),
    ("VHL ligand scaffold",       ["OC(C(=O))NC"],                                            8),
    # Heterocycles of interest
    ("Pyrazolo[3,4-d]pyrimidine", ["c1ncc2[nH]ncc2n1"],                                     8),
    ("Imidazo[1,2-a]pyridine",   ["c1ccc2nccn2c1"],                                          7),
    ("Azaindole",                 ["c1ccc2[nH]cnc2c1"],                                       7),
    ("Indazole",                  ["c1ccc2[nH]ncc2c1"],                                       7),
    ("Macrocycle",                ["[r14,r15,r16,r17,r18,r19,r20]"],                          9),
    # Boronic acids
    ("Boronic acid warhead",      ["[#6]B(O)O"],                                              9),
    # Phosphorus
    ("Phosphonamide",             ["P(=O)([NX3])"],                                            7),
    ("Phosphonate",               ["P(=O)(O)O"],                                               7),
    # Deuterium
    ("Deuterium label",           ["[2H]"],                                                    8),
    # Unnatural amino acids
    ("Unnatural amino acid",      ["[NX3][CX4][CX3](=O)O"],                                  7),
]

# Threshold: require ≥ this many distinct patents before surfacing
MIN_PHARMA_SUPPORT = 2

# Ubiquitous fragments to skip (too common to be informative)
UBIQUITOUS_SMARTS = [
    "c1ccccc1",    # benzene
    "c1ccncc1",    # pyridine
    "CC(=O)",      # acetyl
    "C(=O)N",      # amide (generic)
    "c1ccco1",     # furan
    "c1ccc[nH]1",  # pyrrole
    "C(F)(F)F",    # CF3 (very common)
    "OC(=O)",      # ester/acid
    "S(=O)(=O)",   # sulfone (generic)
]

TOP_N_SIGNALS = 20


def _rdkit_available() -> bool:
    try:
        from rdkit import Chem  # noqa
        return True
    except ImportError:
        return False


def _smarts_match(mol, smarts_list: list[str]) -> bool:
    """Return True if mol matches any SMARTS in the list."""
    try:
        from rdkit import Chem
        for sma in smarts_list:
            patt = Chem.MolFromSmarts(sma)
            if patt and mol.HasSubstructMatch(patt):
                return True
    except Exception:
        pass
    return False


def _is_ubiquitous(mol) -> bool:
    return _smarts_match(mol, UBIQUITOUS_SMARTS)


def _classify_fragment(smiles: str) -> Optional[str]:
    """Return the highest-priority BB class matching this fragment, or None."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if _is_ubiquitous(mol):
            return None
        best_name   = None
        best_pri    = -1
        for name, smarts_list, priority in BB_SMARTS:
            if priority > best_pri and _smarts_match(mol, smarts_list):
                best_name = name
                best_pri  = priority
        return best_name
    except Exception:
        return None


def _brics_fragments(smiles: str) -> list[str]:
    """BRICS decomposition → list of fragment canonical SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import BRICS
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        frags = BRICS.BRICSDecompose(mol)
        cleaned = []
        for f in frags:
            # Remove [*] dummy atoms
            f_clean = re.sub(r'\[\d*\*\]', '', f)
            m = Chem.MolFromSmiles(f_clean)
            if m and m.GetNumHeavyAtoms() >= 4:
                cleaned.append(Chem.MolToSmiles(m))
        return list(set(cleaned))
    except Exception:
        return []


def _make_headline(bb_class: str, assignees: list[str], targets: list[str]) -> str:
    """
    Build a human-readable headline like:
      "Fluorinated cyclopropane in J&J filing targeting KRAS"
    """
    parts = [bb_class]
    if assignees:
        co = assignees[0]
        parts.append(f"in {co} filing")
    if targets:
        tgt = targets[0]
        parts.append(f"targeting {tgt}")
    return " ".join(parts)


def detect_novel_blocks(
    pharma_records: list[dict],
    background_records: list[dict] | None = None,
) -> list[dict]:
    """
    Main entry point.

    Args:
        pharma_records:    Records from sources/pharma_patents.py (source='patent',
                           each has assignee + targets).
        background_records: Full corpus records used to estimate background
                           fragment frequency (for novelty score). Optional.

    Returns:
        pharma_signals list, sorted by support desc, capped at TOP_N_SIGNALS.
    """
    if not _rdkit_available():
        logger.warning("novel_blocks: RDKit not available, skipping.")
        return []

    if not pharma_records:
        return []

    # ── Build background frequency map ──────────────────────────────────────
    bg_freq: dict[str, int] = defaultdict(int)
    if background_records:
        for rec in background_records:
            for frag in _brics_fragments(rec.get("smiles", "")):
                bg_freq[frag] += 1

    # ── Aggregate pharma fragment signals ────────────────────────────────────
    # key: (bb_class, fragment_smiles) → aggregated signal
    agg: dict[tuple[str, str], dict] = {}

    for rec in pharma_records:
        smiles   = rec.get("smiles", "")
        assignee = rec.get("assignee") or "Unknown"
        targets  = rec.get("targets") or []
        date     = rec.get("date") or "2024"
        ref_id   = rec.get("ref_id") or ""

        frags = _brics_fragments(smiles)
        if not frags:
            # Try the full molecule if BRICS yields nothing
            frags = [smiles]

        for frag in frags:
            bb_class = _classify_fragment(frag)
            if not bb_class:
                continue

            key = (bb_class, frag)
            if key not in agg:
                agg[key] = {
                    "bb_class":      bb_class,
                    "fragment_smiles": frag,
                    "assignees":     set(),
                    "targets":       set(),
                    "patent_ids":    set(),
                    "dates":         [],
                    "bg_count":      bg_freq.get(frag, 0),
                }
            sig = agg[key]
            sig["assignees"].add(assignee)
            sig["targets"].update(targets)
            sig["patent_ids"].add(ref_id)
            sig["dates"].append(date)

    # ── Filter, score, sort ──────────────────────────────────────────────────
    signals = []
    for (bb_class, frag), sig in agg.items():
        support = len(sig["patent_ids"])
        if support < MIN_PHARMA_SUPPORT:
            continue

        # Novelty: high if fragment is rare in background corpus
        bg = sig["bg_count"]
        novelty_score = 1.0 / (1.0 + bg)

        # Combined score: support × novelty
        score = round(support * novelty_score * 10, 1)

        assignees_list = sorted(sig["assignees"])
        targets_list   = sorted(sig["targets"])
        dates          = sorted(sig["dates"])
        first_seen     = dates[0] if dates else "2024"

        # Sparkline over years (simple counts per year, last 5)
        year_counts: dict[str, int] = defaultdict(int)
        for d in sig["dates"]:
            year_counts[d[:4]] += 1
        all_years = sorted(year_counts.keys())[-5:]
        sparkline = [year_counts.get(y, 0) for y in all_years]

        signals.append({
            "bb_class":        bb_class,
            "fragment_smiles": frag,
            "assignees":       assignees_list,
            "targets":         targets_list,
            "support":         support,
            "score":           score,
            "novelty":         round(novelty_score, 3),
            "bg_count":        bg,
            "headline":        _make_headline(bb_class, assignees_list, targets_list),
            "badge":           "novel" if novelty_score > 0.5 else "emerging",
            "first_seen":      first_seen,
            "sparkline":       sparkline,
            "years":           all_years,
        })

    # Sort by score desc, then support desc
    signals.sort(key=lambda s: (-s["score"], -s["support"]))
    signals = signals[:TOP_N_SIGNALS]

    # Add sequential IDs
    for i, sig in enumerate(signals, 1):
        sig["id"] = f"ps{i}"

    logger.info("novel_blocks: %d pharma signals detected", len(signals))
    return signals
