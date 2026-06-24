"""
fragment.py — BRICS track: discovery feed for unnamed specific blocks.

Decomposes each molecule into building-block-like synthons via RDKit BRICS,
strips dummy atoms, canonicalizes, counts molecules-containing-block (set, not
per-cut). This finds *specific* emergent blocks that have no named class yet.

The CLASS track (classify.py) is primary for the feed. BRICS is a secondary
discovery layer surfacing un-named blocks that are rising.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Optional

try:
    from rdkit import Chem
    from rdkit.Chem.BRICS import BRICSDecompose
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False


MIN_HEAVY_ATOMS = 4


def brics_fragments(smiles: str) -> set[str]:
    """
    Return the set of canonical BRICS fragment SMILES for a molecule.
    Strips dummy [*] atoms. Skips fragments with < MIN_HEAVY_ATOMS heavy atoms.
    Returns empty set on failure.
    """
    if not RDKIT_OK or not smiles:
        return set()
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return set()
        raw = BRICSDecompose(mol, minFragmentSize=MIN_HEAVY_ATOMS)
        result = set()
        for f in raw:
            # strip dummy atoms [*] and re-canonicalize
            stripped = _strip_dummies(f)
            if stripped:
                result.add(stripped)
        return result
    except Exception:
        return set()


def _strip_dummies(smiles: str) -> Optional[str]:
    """Remove [*] dummy atoms from a BRICS fragment, return canonical SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Remove dummy atoms (atomic num 0)
        em = Chem.RWMol(mol)
        atoms_to_remove = [a.GetIdx() for a in em.GetAtoms() if a.GetAtomicNum() == 0]
        for idx in sorted(atoms_to_remove, reverse=True):
            em.RemoveAtom(idx)
        try:
            Chem.SanitizeMol(em)
        except Exception:
            return None
        # Must have >= MIN_HEAVY_ATOMS heavy atoms after stripping
        if em.GetNumHeavyAtoms() < MIN_HEAVY_ATOMS:
            return None
        canon = Chem.MolToSmiles(em, canonical=True)
        return canon if canon else None
    except Exception:
        return None


def build_fragment_groups(records: list[dict]) -> dict[str, dict]:
    """
    Given standardized records (must have 'smiles', 'source', 'date', 'assignee'),
    decompose every molecule and bucket each BRICS fragment.

    Returns a dict keyed by fragment canonical SMILES:
      {
        frag_smiles: {
          "smiles": str,
          "members": set of (inchikey, source) tuples,
          "dates": list of str,
          "sources": set of str,
          "assignees": set of str,
        }
      }
    """
    groups: dict[str, dict] = {}

    for rec in records:
        smi = rec.get("smiles", "")
        frags = brics_fragments(smi)
        for f in frags:
            if f not in groups:
                groups[f] = {
                    "smiles": f,
                    "members": set(),
                    "dates": [],
                    "sources": set(),
                    "assignees": set(),
                }
            g = groups[f]
            ik = rec.get("inchikey", smi)
            src = rec.get("source", "unknown")
            g["members"].add((ik, src))
            if rec.get("date"):
                g["dates"].append(rec["date"])
            g["sources"].add(src)
            if rec.get("assignee"):
                g["assignees"].add(rec["assignee"])

    return groups
