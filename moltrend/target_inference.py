"""
target_inference.py — Extract drug targets from patent title + abstract text.

Two modes:
  1. Keyword dictionary (always runs) — covers ~95% of common oncology/immunology targets.
     Returns all target names that appear in the text.
  2. LLM fallback (opt-in, ANTHROPIC_API_KEY) — for ambiguous text or novel programs.
     Sends title + abstract to Haiku, asks for structured target list.

Usage:
    from moltrend.target_inference import infer_targets
    targets = infer_targets(title="...", abstract="...")
    # → ["KRAS G12C", "SOS1"]
"""
from __future__ import annotations

import json
import logging
import os
import re
from typing import Optional

logger = logging.getLogger(__name__)


# ── Keyword target dictionary ─────────────────────────────────────────────────
# Maps regex patterns → canonical target name.
# Order matters: more specific patterns first.

TARGET_PATTERNS: list[tuple[str, str]] = [
    # Oncology — kinases
    (r"\bKRAS\s*G12[CDVRS]\b",         "KRAS G12C/D/V"),
    (r"\bKRAS\b",                       "KRAS"),
    (r"\bNRAS\b",                       "NRAS"),
    (r"\bHRAS\b",                       "HRAS"),
    (r"\bBRAF\s*V600",                  "BRAF V600E"),
    (r"\bBRAF\b",                       "BRAF"),
    (r"\bEGFR\s*(C797S|exon\s*20|L858R)?", "EGFR"),
    (r"\bHER2\b|\bERBB2\b",            "HER2"),
    (r"\bHER3\b|\bERBB3\b",            "HER3"),
    (r"\bALK\b",                        "ALK"),
    (r"\bRET\b",                        "RET"),
    (r"\bROS1\b",                       "ROS1"),
    (r"\bFGFR[1-4]?\b",                "FGFR"),
    (r"\bMET\b|\bc-Met\b",             "MET"),
    (r"\bVEGFR[1-3]?\b",              "VEGFR"),
    (r"\bPDGFR[AB]?\b",               "PDGFR"),
    (r"\bSRC\b",                        "SRC"),
    (r"\bABL1?\b",                      "ABL"),
    (r"\bBCR.?ABL\b",                  "BCR-ABL"),
    (r"\bBTK\b",                        "BTK"),
    (r"\bITK\b",                        "ITK"),
    (r"\bTEC\b",                        "TEC"),
    (r"\bPI3K[αβγδ]?\b|\bPI3K\b",    "PI3K"),
    (r"\bAKT[123]?\b|\bPKB\b",        "AKT"),
    (r"\bmTOR\b|\bTORC[12]\b",         "mTOR"),
    (r"\bMEK[12]?\b|\bMAP2K[12]\b",  "MEK"),
    (r"\bERK[12]?\b|\bMAPK[13]\b",   "ERK"),
    (r"\bCDK[1-9](?:\/[1-9])?\b",    "CDK"),
    (r"\bCDK4.?6\b",                   "CDK4/6"),
    (r"\bAURK[AB]?\b|\bAurora\b",     "Aurora kinase"),
    (r"\bPLK1?\b",                     "PLK1"),
    (r"\bWEE1\b",                      "WEE1"),
    (r"\bCHK[12]\b",                   "CHK1/2"),
    (r"\bATR\b",                       "ATR"),
    (r"\bATM\b",                       "ATM"),
    (r"\bPARP[12]?\b",                 "PARP"),
    (r"\bSOS1?\b",                     "SOS1"),
    (r"\bSHP2\b|\bPTPN11\b",          "SHP2"),
    (r"\bFAK\b|\bPTK2\b",             "FAK"),
    (r"\bLCK\b",                       "LCK"),
    (r"\bZAP70\b",                     "ZAP-70"),
    (r"\bIRAK[14]?\b",                "IRAK4"),
    (r"\bRIP[K]?[123]\b|\bRIPK\b",   "RIPK"),
    (r"\bTBK1\b",                      "TBK1"),
    (r"\bIKK[βε]?\b|\bIKBK\b",       "IKK"),
    (r"\bSTING\b|\bTMEM173\b",        "STING"),
    (r"\bJAK[1-3]?\b",                "JAK"),
    (r"\bTYK2\b",                      "TYK2"),
    (r"\bSTAT[1-6]\b",                "STAT"),
    # Oncology — E3 ligases / degraders / CELMoDs
    (r"\bcereblon\b|\bCRBN\b|\bIMiD\b|\bCELMoD\b", "CRBN (cereblon)"),
    (r"\bVHL\b|\bvon Hippel.?Lindau\b", "VHL"),
    (r"\bMDM2\b",                      "MDM2"),
    (r"\bIAP\b|\bXIAP\b|\bBIRC\b",   "IAP"),
    (r"\bPROTAC\b|\bTPD\b|\bdegrader\b|\btargeted protein degradation\b", "TPD (degrader)"),
    (r"\bglue\b.{0,20}molecul|\bmolecular glue\b", "Molecular glue"),
    (r"\bIKZF[123]\b|\bAiolos\b|\bIkaros\b", "IKZF1/3"),
    (r"\bGSPT1\b|\beIF3e\b",          "GSPT1"),
    (r"\bCK1[αβ]?\b",                 "CK1α"),
    # Oncology — transcription factors / PPI
    (r"\bBRD[234]?\b|\bBET\b|\bbromodomain\b", "BRD4 / BET"),
    (r"\bMYC\b|\bc-Myc\b",            "MYC"),
    (r"\bBCL.?2\b",                   "BCL-2"),
    (r"\bBCL.?XL\b|\bBCL2L1\b",     "BCL-XL"),
    (r"\bMCL.?1\b",                   "MCL-1"),
    (r"\bp53\b|\bTP53\b",             "p53"),
    (r"\bMDM4\b|\bMDMX\b",           "MDM4"),
    (r"\bWNT\b|\bβ.?catenin\b|\bCTNNB1\b", "β-catenin / WNT"),
    (r"\bHedgehog\b|\bSMO\b|\bGLI[123]?\b", "Hedgehog / SMO"),
    (r"\bNotch\b",                     "Notch"),
    # Epigenetics
    (r"\bEZH2\b",                      "EZH2"),
    (r"\bDOT1L\b",                    "DOT1L"),
    (r"\bLSD1\b|\bKDM1A\b",          "LSD1"),
    (r"\bHDAC[1-9]?\b",              "HDAC"),
    (r"\bBET\b",                       "BET bromodomain"),
    (r"\bPRMT[157]\b",               "PRMT"),
    # Immunology / inflammation
    (r"\bIL.?[0-9]+\b|\binterleukin\b", "IL receptor"),
    (r"\bTNF[αR]?\b|\btumour necrosis\b", "TNFα"),
    (r"\bIL.?17[AR]?\b",             "IL-17"),
    (r"\bIL.?23\b|\bIL.?12\b",      "IL-12/23"),
    (r"\bIFN[αβγ]?\b|\binterferon\b", "IFN"),
    (r"\bPD.?[L1]\b|\bPD-1\b|\bPDCD1\b", "PD-1 / PD-L1"),
    (r"\bCTLA.?4\b",                  "CTLA-4"),
    (r"\bLAG.?3\b",                   "LAG-3"),
    (r"\bTIM.?3\b",                   "TIM-3"),
    (r"\bTIGIT\b",                    "TIGIT"),
    (r"\bCD47\b",                     "CD47"),
    (r"\bCD20\b",                     "CD20"),
    (r"\bCD19\b",                     "CD19"),
    (r"\bBCMA\b|\bTNFRSF17\b",       "BCMA"),
    (r"\bGPCR\b|\bG.protein.coupled\b", "GPCR"),
    (r"\bAdenosine\b|\bA2A[R]?\b|\bA2B[R]?\b", "Adenosine receptor"),
    (r"\bGLP.?1[R]?\b|\bsemaglutide\b|\bliraglutide\b", "GLP-1R"),
    (r"\bGIP[R]?\b",                  "GIPR"),
    # Ion channels / CNS
    (r"\bNav1\.[1-9]\b|\bsodium channel\b", "Nav channel"),
    (r"\bKv[0-9]\b|\bpotassium channel\b", "Kv channel"),
    (r"\bNMDA\b|\bAMPA\b|\bGluR\b",  "Glutamate receptor"),
    (r"\bGABA[AB]?\b",               "GABA receptor"),
    (r"\bdopamine\b|\bDRD[1-5]\b",   "Dopamine receptor"),
    (r"\bserotonin\b|\b5.HT[1-7]\b|\bHTR\b", "Serotonin receptor"),
    # Proteases / other enzymes
    (r"\bDPP.?4\b|\bCD26\b",         "DPP-4"),
    (r"\bACE2?\b",                    "ACE"),
    (r"\bcaspase.?[1-9]\b",          "Caspase"),
    (r"\bMatriptase\b|\bST14\b",     "Matriptase"),
    # Antiviral
    (r"\bSARS.?CoV.?2\b|\bCOVID\b|\bnSP[0-9]+\b|\bMpro\b|\b3CL\b", "SARS-CoV-2"),
    (r"\bHIV\b|\bprotease\b.{0,10}\bHIV\b|\bintegrase\b", "HIV"),
    (r"\bHCV\b|\bhepatitis C\b|\bNS[35]B?\b",            "HCV"),
    (r"\bRSV\b|\brespiratory syncytial\b",                "RSV"),
    (r"\bInfluenza\b|\bNeuraminidase\b|\bHA\b.{0,10}influenza", "Influenza"),
    # Antibacterial
    (r"\bDNA gyrase\b|\btopoisomerase\b|\bGyrB\b",        "DNA Gyrase"),
    (r"\bFabI\b|\benoyl.ACP\b",      "FabI"),
]

# Compile once at import
_COMPILED_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(pat, re.IGNORECASE), name)
    for pat, name in TARGET_PATTERNS
]

# Common ubiquitous terms to exclude even if they match
_EXCLUDE = {"a", "an", "the", "in", "of", "and", "or", "for", "to", "as"}


def infer_targets(title: str = "", abstract: str = "", use_llm: bool = True) -> list[str]:
    """
    Infer drug targets from patent title + abstract text.
    Returns a deduplicated list of canonical target names.
    """
    text = f"{title} {abstract}"
    if not text.strip():
        return []

    # ── Keyword pass ──────────────────────────────────────────────────────────
    found: dict[str, bool] = {}
    for pattern, name in _COMPILED_PATTERNS:
        if pattern.search(text):
            found[name] = True

    targets = list(found.keys())

    # ── LLM pass (if keyword pass returns nothing useful and key available) ───
    if not targets and use_llm and os.environ.get("ANTHROPIC_API_KEY"):
        targets = _llm_targets(title, abstract[:800])

    return targets[:6]  # cap at 6 per patent


def _llm_targets(title: str, abstract_snippet: str) -> list[str]:
    """Ask Haiku for targets when keyword pass fails. Returns list of strings."""
    try:
        import anthropic
        client = anthropic.Anthropic()
        msg = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=120,
            system=(
                "You extract drug targets from patent text. "
                "Return ONLY a JSON array of target name strings (e.g. [\"KRAS\",\"SOS1\"]). "
                "If no clear target, return []. No explanation."
            ),
            messages=[{"role": "user", "content": f"Title: {title}\n\nAbstract: {abstract_snippet}"}],
        )
        text = msg.content[0].text.strip()
        start = text.find("[")
        end   = text.rfind("]") + 1
        if start >= 0 and end > start:
            result = json.loads(text[start:end])
            return [str(t) for t in result if t]
    except Exception as e:
        logger.debug("LLM target inference failed: %s", e)
    return []


def infer_modality(title: str = "", abstract: str = "") -> Optional[str]:
    """
    Infer drug modality from text.
    Returns one of: 'small_molecule', 'antibody', 'adc', 'protac', 'peptide', 'other', None
    """
    text = f"{title} {abstract}".lower()
    if any(k in text for k in ["protac", "degrader", "molecular glue", "celmod"]):
        return "protac_glue"
    if any(k in text for k in ["antibody-drug conjugate", "adc", " adc ", "linker-payload"]):
        return "adc"
    if any(k in text for k in ["antibody", "mab", "monoclonal", "bispecific"]):
        return "antibody"
    if any(k in text for k in ["peptide", "cyclic peptide", "macrolide", "stapled"]):
        return "peptide"
    if any(k in text for k in ["oligonucleotide", "sirna", "mrna", "antisense"]):
        return "nucleic_acid"
    return "small_molecule"
