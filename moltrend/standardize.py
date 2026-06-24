"""
standardize.py — Single trusted gate for all SMILES entering MolTrend.

Every structure from every source passes through here:
  parse → sanitize → desalt (largest fragment) → neutralize →
  canonical SMILES + InChIKey + descriptors

The InChIKey is the dedup key. Dedup is per (inchikey, source) so a molecule
in both a patent and an approval counts once per source — that split is signal.
"""
from __future__ import annotations

from typing import Optional
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, inchi
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False
    logging.warning("RDKit not available — standardize will pass through raw SMILES")

logger = logging.getLogger(__name__)


def standardize(smiles: str) -> Optional[dict]:
    """
    Standardize a SMILES string. Returns a dict on success, None on failure.

    Output dict keys:
      smiles      canonical SMILES
      inchikey    InChIKey (dedup key)
      mw          molecular weight (Da)
      clogp       Crippen cLogP
      hbd         H-bond donors
      hba         H-bond acceptors
      tpsa        topological polar surface area
      rotb        rotatable bonds
    """
    if not smiles or not isinstance(smiles, str):
        return None
    smiles = smiles.strip()
    if not smiles:
        return None

    if not RDKIT_OK:
        # Passthrough mode — no descriptors
        return {
            "smiles": smiles,
            "inchikey": _fake_inchikey(smiles),
            "mw": 0.0, "clogp": 0.0, "hbd": 0, "hba": 0, "tpsa": 0.0, "rotb": 0,
        }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Desalt: keep largest fragment
        remover = rdMolStandardize.FragmentRemover()
        mol = remover.remove(mol)
        if mol is None:
            return None

        # Sanitize
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return None

        # Neutralize common charges
        try:
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)
        except Exception:
            pass  # non-critical

        # Canonical SMILES (round-trip validation)
        canon = Chem.MolToSmiles(mol, canonical=True)
        if not canon:
            return None
        mol2 = Chem.MolFromSmiles(canon)
        if mol2 is None:
            return None

        # InChIKey
        try:
            ik = inchi.MolToInchiKey(mol2)
        except Exception:
            ik = None
        if not ik:
            return None

        # Descriptors
        mw   = round(Descriptors.MolWt(mol2), 2)
        clogp = round(Descriptors.MolLogP(mol2), 2)
        hbd  = rdMolDescriptors.CalcNumHBD(mol2)
        hba  = rdMolDescriptors.CalcNumHBA(mol2)
        tpsa = round(rdMolDescriptors.CalcTPSA(mol2), 1)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol2)

        return {
            "smiles": canon,
            "inchikey": ik,
            "mw": mw,
            "clogp": clogp,
            "hbd": hbd,
            "hba": hba,
            "tpsa": tpsa,
            "rotb": rotb,
        }

    except Exception as e:
        logger.debug("standardize failed for %s: %s", smiles[:50], e)
        return None


def morgan_fp(smiles: str, radius: int = 2, nbits: int = 1024) -> Optional[list]:
    """Return Morgan fingerprint as a list of ints, or None."""
    if not RDKIT_OK:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        return list(fp)
    except Exception:
        return None


def _fake_inchikey(smiles: str) -> str:
    """Deterministic placeholder InChIKey when RDKit unavailable."""
    import hashlib
    h = hashlib.sha256(smiles.encode()).hexdigest()[:27].upper()
    return f"{h[:14]}-{h[14:24]}-N"
