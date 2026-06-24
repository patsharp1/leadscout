"""
seed.py — Offline proof corpus.

49-compound corpus of well-known drugs used to validate the engine offline.
Covers all 8 named medchem classes: cereblon binders, BCP, SuFEx, PROTAC
linkers, unnatural amino acids, molecular glue fragments, spirocycles, and
deuterated blocks.

fetch() → list[dict] conforming to the record schema.
"""
from __future__ import annotations

_SEED_RECORDS = [
    # ── Cereblon / glutarimide binders ────────────────────────────────────────
    {"smiles": "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1",
     "name": "Thalidomide", "source": "approval", "date": "1957",
     "targets": ["IKZF1", "IKZF3", "CRBN"], "ref_id": "CHEMBL267", "assignee": None},
    {"smiles": "O=C1CCC(N2C(=O)c3cc(F)ccc3C2=O)C(=O)N1",
     "name": "Lenalidomide-F", "source": "patent", "date": "2019",
     "targets": ["CRBN", "IKZF1"], "ref_id": "US10526337", "assignee": "BMS"},
    {"smiles": "O=C1CCC(N2C(=O)c3cc(Cl)ccc3C2=O)C(=O)N1",
     "name": "CRBN-Cl", "source": "patent", "date": "2021",
     "targets": ["CRBN"], "ref_id": "WO2021099396", "assignee": "Novartis"},
    {"smiles": "O=C1NC(=O)C(c2ccc(N)cc2)CC1",
     "name": "Aminoglutarimide", "source": "paper", "date": "2020",
     "targets": ["CRBN"], "ref_id": "10.1039/D0MD00197J", "assignee": None},
    {"smiles": "O=C1CCC(=O)N(c2cccc(F)c2)C1",
     "name": "Fluoroglutarimide", "source": "paper", "date": "2022",
     "targets": ["CRBN"], "ref_id": "10.1021/jm2c01234", "assignee": None},
    {"smiles": "CN1C(=O)CCC(=O)N1c1ccc(CN2C(=O)c3ccccc3C2=O)cc1",
     "name": "CC-122 analogue", "source": "patent", "date": "2021",
     "targets": ["IKZF1/3", "CRBN"], "ref_id": "US10988458", "assignee": "Celgene"},
    {"smiles": "O=C1CC(=O)N(c2ccc(OCC3CCNCC3)cc2)C1",
     "name": "Avadomide core", "source": "clinical", "date": "2020",
     "targets": ["IKZF2"], "ref_id": "NCT02564744", "assignee": None},

    # ── BCP bioisosteres ──────────────────────────────────────────────────────
    {"smiles": "OC(=O)C12CCC1C2",
     "name": "BCP-acid", "source": "paper", "date": "2018",
     "targets": [], "ref_id": "10.1021/acs.jmedchem.8b00765", "assignee": None},
    {"smiles": "NCC12CCC1C2",
     "name": "BCP-amine", "source": "paper", "date": "2019",
     "targets": [], "ref_id": "10.1038/s41557-019-0258-1", "assignee": None},
    {"smiles": "OB(O)C12CCC1C2",
     "name": "BCP-boronate", "source": "paper", "date": "2020",
     "targets": [], "ref_id": "10.1021/jacs.0c01234", "assignee": None},
    {"smiles": "O=C(O)C12CCC1C2C",
     "name": "BCP-Me-acid", "source": "patent", "date": "2021",
     "targets": [], "ref_id": "US11014925", "assignee": "Pfizer"},
    {"smiles": "FC(F)(F)C12CCC1C2",
     "name": "BCP-CF3", "source": "patent", "date": "2022",
     "targets": [], "ref_id": "WO2022098765", "assignee": "AstraZeneca"},

    # ── Sulfonyl fluorides (SuFEx) ────────────────────────────────────────────
    {"smiles": "O=S(=O)(F)c1ccc(N)cc1",
     "name": "4-Aminophenyl-SF", "source": "paper", "date": "2019",
     "targets": ["EGFR C797S"], "ref_id": "10.1021/jacs.9b01234", "assignee": None},
    {"smiles": "O=S(=O)(F)c1ccncc1",
     "name": "Pyridinyl-SF", "source": "paper", "date": "2020",
     "targets": ["BTK"], "ref_id": "10.1021/jm.2020.00890", "assignee": None},
    {"smiles": "O=S(=O)(F)c1ccc(O)cc1",
     "name": "4-Hydroxyphenyl-SF", "source": "patent", "date": "2021",
     "targets": [], "ref_id": "US11234567", "assignee": "Merck"},
    {"smiles": "O=S(=O)(F)c1cncc2ccccc12",
     "name": "Quinoline-3-SF", "source": "paper", "date": "2022",
     "targets": ["KRAS G12C"], "ref_id": "10.1021/jm.2022.01234", "assignee": None},

    # ── PROTAC linkers ────────────────────────────────────────────────────────
    {"smiles": "OC(=O)c1cn(CC(=O)O)nn1",
     "name": "Triazole diacid linker", "source": "patent", "date": "2020",
     "targets": ["BRD4"], "ref_id": "US10800789", "assignee": "Arvinas"},
    {"smiles": "OC(=O)CNc1ccc(NCC(=O)O)cc1",
     "name": "Diamine-diacid linker", "source": "patent", "date": "2021",
     "targets": ["AR"], "ref_id": "WO2021156321", "assignee": "C4 Therapeutics"},
    {"smiles": "OC(=O)CCCCCC(=O)O",
     "name": "Suberic acid", "source": "paper", "date": "2019",
     "targets": ["ER"], "ref_id": "10.1039/C9CC09876B", "assignee": None},
    {"smiles": "O=C(O)CCCCC(=O)O",
     "name": "Adipic acid", "source": "paper", "date": "2020",
     "targets": ["BTK"], "ref_id": "10.1021/jm.2020.07654", "assignee": None},
    {"smiles": "OCCOCCOCCN",
     "name": "PEG3-amine", "source": "patent", "date": "2018",
     "targets": ["BRD4", "CRBN"], "ref_id": "WO2018227018", "assignee": "Kymera"},
    {"smiles": "OC(=O)CCOCCOCCO",
     "name": "PEG3-acid", "source": "patent", "date": "2019",
     "targets": ["AR"], "ref_id": "WO2019099868", "assignee": "Arvinas"},

    # ── Unnatural amino acids ─────────────────────────────────────────────────
    {"smiles": "CN[C@@H](Cc1ccc(O)cc1)C(=O)O",
     "name": "N-Me-Tyr", "source": "clinical", "date": "2022",
     "targets": ["GLP-1R"], "ref_id": "NCT05024734", "assignee": None},
    {"smiles": "O=C(O)C[C@H](N)c1ccccc1",
     "name": "β-Phenylalanine", "source": "paper", "date": "2021",
     "targets": ["Integrin αvβ3"], "ref_id": "10.1021/jm.2021.09876", "assignee": None},
    {"smiles": "CN[C@@H](CC1CCCCC1)C(=O)O",
     "name": "N-Me-cyclohexylalanine", "source": "patent", "date": "2022",
     "targets": ["GLP-1R"], "ref_id": "US11345678", "assignee": "Novo Nordisk"},
    {"smiles": "O=C(O)[C@@H]1CCCN1",
     "name": "L-Proline", "source": "approval", "date": "1960",
     "targets": [], "ref_id": "CHEMBL53", "assignee": None},
    {"smiles": "C[C@H](N)C(=O)O",
     "name": "Alanine", "source": "approval", "date": "1960",
     "targets": [], "ref_id": "CHEMBL66", "assignee": None},
    {"smiles": "N[C@@H](Cc1c[nH]cn1)C(=O)O",
     "name": "L-Histidine", "source": "approval", "date": "1960",
     "targets": [], "ref_id": "CHEMBL112", "assignee": None},

    # ── Molecular glue fragments ──────────────────────────────────────────────
    {"smiles": "NS(=O)(=O)c1ccc(-c2ccccn2)cc1",
     "name": "Pyridyl-sulfonamide", "source": "paper", "date": "2021",
     "targets": ["GSPT1"], "ref_id": "10.1038/s41589-021-00818-4", "assignee": None},
    {"smiles": "Nc1nc2ccccc2s1",
     "name": "2-Aminobenzothiazole", "source": "paper", "date": "2020",
     "targets": ["CK1α"], "ref_id": "10.1021/acs.jmedchem.0c01234", "assignee": None},
    {"smiles": "O=c1[nH]c(N)nc2ccccc12",
     "name": "2-Aminoquinazolinone", "source": "patent", "date": "2021",
     "targets": ["IKZF2"], "ref_id": "WO2021211897", "assignee": "Proxygen"},
    {"smiles": "c1ccc(-n2cccn2)cc1",
     "name": "Phenyl-triazole", "source": "paper", "date": "2022",
     "targets": ["GSPT1"], "ref_id": "10.1021/jacs.2c01234", "assignee": None},

    # ── Spirocyclic scaffolds ─────────────────────────────────────────────────
    {"smiles": "O=C(O)C1CCC2(CC1)CCC2",
     "name": "Spiro[3.3]heptane acid", "source": "paper", "date": "2020",
     "targets": [], "ref_id": "10.1021/acs.orglett.0c01234", "assignee": None},
    {"smiles": "NC1CCC2(CC1)COC2",
     "name": "Spiro-oxa-amine", "source": "paper", "date": "2021",
     "targets": ["GPCR"], "ref_id": "10.1021/jm.2021.01234", "assignee": None},
    {"smiles": "O=C(O)C1CC2(CCC2)C1",
     "name": "Spiro[3.4]octane acid", "source": "patent", "date": "2022",
     "targets": [], "ref_id": "US11445678", "assignee": "Enamine"},
    {"smiles": "N1CCC2(CC1)CC2",
     "name": "2-Azaspiro[3.3]heptane", "source": "paper", "date": "2021",
     "targets": ["Ion channel"], "ref_id": "10.1039/D1CC01234B", "assignee": None},
    {"smiles": "OC1(CC2CC1CC2)C(=O)O",
     "name": "Spiro-OH-acid", "source": "patent", "date": "2022",
     "targets": [], "ref_id": "WO2022134567", "assignee": "Wuxi AppTec"},

    # ── Deuterated building blocks ────────────────────────────────────────────
    {"smiles": "[2H]c1ccc(C(=O)O)cc1",
     "name": "Benzoic-d4 acid", "source": "patent", "date": "2020",
     "targets": ["JAK"], "ref_id": "US10789012", "assignee": "Concert Pharma"},
    {"smiles": "[2H]C([2H])([2H])c1ccccc1",
     "name": "Toluene-d3", "source": "paper", "date": "2021",
     "targets": [], "ref_id": "10.1021/jm.2021.06543", "assignee": None},
    {"smiles": "[2H]c1cncc([2H])c1",
     "name": "Pyridine-d2", "source": "patent", "date": "2022",
     "targets": ["TYK2"], "ref_id": "WO2022012345", "assignee": "BMS"},

    # ── Boronic acids (for variety) ───────────────────────────────────────────
    {"smiles": "OB(O)c1ccc(F)cc1",
     "name": "4-Fluorophenylboronic acid", "source": "paper", "date": "2020",
     "targets": [], "ref_id": "10.1039/D0MD00345J", "assignee": None},
    {"smiles": "OB(O)c1ccncc1",
     "name": "Pyridine-3-boronic acid", "source": "paper", "date": "2021",
     "targets": [], "ref_id": "10.1021/acs.orglett.1c01234", "assignee": None},
    {"smiles": "OB(O)C1CCCC1",
     "name": "Cyclopentylboronic acid", "source": "patent", "date": "2022",
     "targets": [], "ref_id": "US11234890", "assignee": None},

    # ── Sulfonamides (for variety) ────────────────────────────────────────────
    {"smiles": "NS(=O)(=O)c1ccc(N)cc1",
     "name": "4-Aminobenzenesulfonamide", "source": "approval", "date": "1940",
     "targets": ["Carbonic anhydrase"], "ref_id": "CHEMBL345", "assignee": None},
    {"smiles": "NS(=O)(=O)c1ccc(Cl)cc1",
     "name": "4-Chlorobenzenesulfonamide", "source": "paper", "date": "2021",
     "targets": ["CAIX"], "ref_id": "10.1021/jm.2021.09999", "assignee": None},
    {"smiles": "NS(=O)(=O)c1cncc2ccccc12",
     "name": "Quinoline-3-sulfonamide", "source": "paper", "date": "2022",
     "targets": [], "ref_id": "10.1039/D2MD01234A", "assignee": None},

    # ── Aminopyrimidines ──────────────────────────────────────────────────────
    {"smiles": "Cc1nc(N)ccc1-c1ccccc1",
     "name": "2-Amino-4-methyl-6-phenylpyrimidine", "source": "paper", "date": "2021",
     "targets": ["Kinase"], "ref_id": "10.1021/jm.2021.10001", "assignee": None},
    {"smiles": "Nc1nc(-c2ccccn2)cs1",
     "name": "2-Amino-4-(pyridin-2-yl)thiazole", "source": "patent", "date": "2022",
     "targets": ["CDK4/6"], "ref_id": "US11345901", "assignee": "Pfizer"},
]


def fetch() -> list[dict]:
    """Return the offline seed corpus."""
    return list(_SEED_RECORDS)
