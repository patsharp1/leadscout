"""
classify.py — CLASS track (primary): group molecules into named medchem block
classes via SMARTS + spiro detector.

A molecule joins every class it matches. This is what the app surfaces; chemists
reason in classes ("BCP", "cereblon binder"), not BRICS fragments.

Extend CLASS_SMARTS as the corpus grows — it IS the taxonomy users see.
"""
from __future__ import annotations

from typing import Optional

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

# ── Named medchem block classes ──────────────────────────────────────────────
# Key:   display name shown in the feed
# Value: list of SMARTS patterns (molecule matches if ANY pattern matches)

CLASS_SMARTS: dict[str, list[str]] = {
    "Cereblon glutarimide":      ["O=C1CCCC(=O)N1", "O=C1NC(=O)CCC1"],
    "Boronic acid warhead":      ["[#6]B(O)O"],
    "Benzoxaborole":             ["B1OCc2ccccc21"],
    "Sulfonyl fluoride (SuFEx)": ["S(=O)(=O)F"],
    "Acrylamide warhead":        ["[CX3]=[CX3]C(=O)[NX3]", "C#CC(=O)[NX3]"],
    "Aminopyrimidine":           ["[NX3]c1ncccn1", "[NX3]c1ncncc1"],
    "Quinazoline":               ["c1ccc2ncncc2c1"],
    "BCP bioisostere":           ["C12CCC1C2"],
    "Sulfonamide":               ["S(=O)(=O)[NX3]"],
    "Piperazine":                ["C1CNCCN1"],
    "Morpholine":                ["O1CCNCC1"],
    "Benzothiazole":             ["c1ccc2scnc2c1"],
    "Trifluoromethyl":           ["[CX4](F)(F)F"],
    "N-aryl urea":               ["[NX3]C(=O)[NX3]c1ccccc1"],
    # Extended classes ─────────────────────────────────────────────
    "PROTAC linker":             [
        "OC(=O)CCOCCO",           # PEG-acid linker motif
        "N1CCN(CCO)CC1",          # piperazine-PEG
        "OC(=O)c1cn(CC(=O)O)nn1", # triazole diacid
        "OCCOCCOCCN",             # PEG-amine
    ],
    "Unnatural amino acid":      [
        "[NX3][CX4][CX3](=O)[OX2H]",   # general alpha-amino acid
        "C[NX3][CX4][CX3](=O)[OX2H]",  # N-methyl amino acid
    ],
    "Fluorosulfate":             ["OS(=O)(=O)F"],
    "Oxetane":                   ["C1COC1"],
    "Azetidine":                 ["C1CNC1"],
    "Deuterium label":           ["[2H]"],
    "Macrolide / lactone":       ["C1OC(=O)CCCCCC1", "C1OC(=O)CCCCCCC1"],
    "Heterobifunctional linker": ["N1CCN(CC1)C(=O)CCCC(=O)O"],
    # ── Modern warheads & strained ring systems ──────────────────────────────
    # Fluorinated cyclopropanes — rising in KRAS, degrader, and CNS programs
    "Fluorinated cyclopropane":  ["FC1CC1", "C(F)(F)C1CC1", "[F]C1(CC1)"],
    "gem-Difluorocyclopropane":  ["FC1(F)CC1"],
    "Fluorinated cyclobutane":   ["FC1CCC1", "C(F)(F)C1CCC1"],
    "gem-Difluorocyclobutane":   ["FC1(F)CCC1"],
    # Bicyclo[1.1.0]butane — ultra-strained sp3 warhead
    "Bicyclo[1.1.0]butane (BCB)": ["C12CC1C2"],
    # Covalent warheads beyond acrylamide
    "Cyanoacrylamide warhead":   ["N#CC=CC(=O)[NX3]"],
    "Chloroacetamide warhead":   ["ClCC(=O)[NX3]"],
    "Vinylsulfonamide warhead":  ["[CX3]=[CX3]S(=O)(=O)[NX3]"],
    "Epoxide warhead":           ["[CX4]1O[CX4]1"],
    # VHL handle (complement to glutarimide on VHL-recruiting PROTACs)
    "VHL ligand (hydroxyproline)": [
        "OC1CCNC1C(=O)",           # trans-4-hydroxyproline motif
        "[C@@H]1(O)CCN[C@H]1C(=O)",
    ],
    # Phosphorous-based fragments
    "Phosphonamide":             ["[NX3]P(=O)"],
    "Phosphonate ester":         ["OP(=O)(OC)OC"],
    # Boron beyond boronic acid
    "MIDA boronate":             ["B1OC(=O)CN(C)CC(=O)O1"],
    "Boronic ester":             ["B1OC(CO)CO1", "B1OCC(C)(C)O1"],
    # Macrocyclic / cyclic peptide mimetics
    "Macrocyclic scaffold":      [
        "[r14]", "[r15]", "[r16]", "[r17]", "[r18]",
    ],
}

# Pre-compile query mols at import time
_COMPILED: dict[str, list] = {}


def _compile_patterns():
    if _COMPILED or not RDKIT_OK:
        return
    for cls_name, patterns in CLASS_SMARTS.items():
        qmols = []
        for pat in patterns:
            try:
                qm = Chem.MolFromSmarts(pat)
                if qm:
                    qmols.append(qm)
            except Exception:
                pass
        _COMPILED[cls_name] = qmols


def _is_spiro(mol) -> bool:
    """Return True if mol contains any two SSSR rings sharing exactly one atom."""
    try:
        ri = mol.GetRingInfo()
        rings = [set(r) for r in ri.AtomRings()]
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                shared = rings[i] & rings[j]
                if len(shared) == 1:
                    return True
    except Exception:
        pass
    return False


def classify_molecule(smiles: str) -> list[str]:
    """
    Return list of class names that match this molecule.
    A molecule can belong to multiple classes.
    """
    if not RDKIT_OK or not smiles:
        return []
    _compile_patterns()
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
    except Exception:
        return []

    matched = []
    for cls_name, qmols in _COMPILED.items():
        for qm in qmols:
            try:
                if mol.HasSubstructMatch(qm):
                    matched.append(cls_name)
                    break  # matched this class — move to next
            except Exception:
                pass

    # Spirocycle detector (code-based, not SMARTS)
    if _is_spiro(mol):
        matched.append("Spirocyclic scaffold")

    return matched


def build_class_groups(records: list[dict]) -> dict[str, dict]:
    """
    Given standardized records, classify every molecule and bucket by class.

    Returns a dict keyed by class name:
      {
        class_name: {
          "bb_class": str,
          "primary_smiles": str,          # first member's SMILES (as representative)
          "members": set of (inchikey, source),
          "member_smiles": list of str,   # up to 5 representative SMILES
          "dates": list of str,
          "sources": set of str,
          "assignees": set of str,
          "targets": set of str,
          "approval_count": int,
          "patent_count": int,
          "paper_count": int,
        }
      }
    """
    groups: dict[str, dict] = {}

    for rec in records:
        smi = rec.get("smiles", "")
        classes = classify_molecule(smi)

        for cls in classes:
            if cls not in groups:
                groups[cls] = {
                    "bb_class": cls,
                    "primary_smiles": smi,
                    "members": set(),
                    "member_smiles": [],
                    "dates": [],
                    "sources": set(),
                    "assignees": set(),
                    "targets": set(),
                    "approval_count": 0,
                    "patent_count": 0,
                    "paper_count": 0,
                }
            g = groups[cls]
            ik = rec.get("inchikey", smi)
            src = rec.get("source", "unknown")
            key = (ik, src)
            if key not in g["members"]:
                g["members"].add(key)
                if len(g["member_smiles"]) < 5:
                    g["member_smiles"].append(smi)
            if rec.get("date"):
                g["dates"].append(rec["date"])
            g["sources"].add(src)
            if rec.get("assignee"):
                g["assignees"].add(rec["assignee"])
            for t in rec.get("targets", []):
                if t:
                    g["targets"].add(t)
            if src == "approval":
                g["approval_count"] += 1
            elif src == "patent":
                g["patent_count"] += 1
            elif src in ("paper", "clinical"):
                g["paper_count"] += 1

    return groups


def list_classes() -> list[str]:
    """Return all defined class names (including spirocycle)."""
    return list(CLASS_SMARTS.keys()) + ["Spirocyclic scaffold"]
