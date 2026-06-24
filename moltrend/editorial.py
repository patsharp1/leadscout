"""
editorial.py — Write card copy for trend cards.

Two paths:
  LLM-enhanced (Haiku): rewrite copy in an editorial "taste" voice using ONLY
    computed numbers — never invents data. Requires ANTHROPIC_API_KEY env var.
  Deterministic fallback: template-based, always runs, no secrets needed.

The fallback is the default; the LLM path is opt-in via env var.
"""
from __future__ import annotations

import os
import json
import logging
from string import Template

logger = logging.getLogger(__name__)

# ── Deterministic templates ───────────────────────────────────────────────────

_ONE_LINER_TEMPLATES = {
    "hot": Template(
        "${cls} appears in ${evidence} structures — accelerating across "
        "${src_str} with ${assignee_note}."
    ),
    "rising": Template(
        "Rising signal: ${cls} in ${evidence} structures, "
        "velocity positive across ${src_str}."
    ),
    "emerging": Template(
        "New entrant: ${cls} first detected in the most recent data window "
        "(${evidence} supporting structures)."
    ),
    "steady": Template(
        "${cls} holds a steady niche — ${evidence} structures, "
        "consistent across ${src_str}."
    ),
}

_DETAIL_TEMPLATES = {
    "hot": Template(
        "${cls} is accelerating. ${evidence} structures across "
        "${src_str} show ${vel_desc} velocity. "
        "${whitespace_note}"
        "${assignee_detail}"
        "Top targets: ${targets_str}."
    ),
    "rising": Template(
        "${cls} shows upward momentum in ${src_str}. "
        "${evidence} structures captured; "
        "${whitespace_note}"
        "${assignee_detail}"
        "Associated targets: ${targets_str}."
    ),
    "emerging": Template(
        "${cls} is a new entrant — first appeared in the latest data window. "
        "${evidence} structures so far. "
        "${whitespace_note}"
        "${assignee_detail}"
    ),
    "steady": Template(
        "${cls} maintains steady presence with ${evidence} structures in "
        "${src_str}. Stable velocity. "
        "${whitespace_note}"
        "${assignee_detail}"
    ),
}

_ACTION_TEMPLATES = {
    "hot":      Template("Stock ${cls} building blocks now — demand is accelerating."),
    "rising":   Template("Monitor ${cls} closely; consider diversifying your ${cls} panel."),
    "emerging": Template("Early signal: build awareness of ${cls} — evaluate for your programs."),
    "steady":   Template("${cls} is a reliable workhorse — maintain standard inventory."),
}


def _src_str(sources: list[str]) -> str:
    mapping = {"patent": "patents", "paper": "literature", "approval": "approved drugs",
               "clinical": "clinical candidates"}
    parts = [mapping.get(s, s) for s in sources]
    if len(parts) == 1:
        return parts[0]
    return ", ".join(parts[:-1]) + " and " + parts[-1]


def _vel_desc(score: float) -> str:
    if score >= 90:
        return "very high"
    if score >= 70:
        return "high"
    if score >= 50:
        return "moderate"
    return "low"


def _whitespace_note(ws: bool, sources: list[str]) -> str:
    if ws and "approval" not in sources:
        return "⚡ Whitespace: patent/literature signal with no approved drugs yet — predictive. "
    return ""


def _assignee_detail(assignees: list[str]) -> str:
    if not assignees:
        return ""
    top = assignees[:3]
    return "Key assignees: " + ", ".join(top) + ". "


def _targets_str(targets: list[str]) -> str:
    if not targets:
        return "not yet assigned"
    return ", ".join(targets[:4])


def write_copy_deterministic(trend: dict) -> dict:
    """Fill one_liner, detail, action using templates. Returns updated trend dict."""
    badge     = trend.get("badge", "steady")
    cls       = trend.get("bb_class", "Unknown")
    evidence  = trend.get("evidence", 0)
    sources   = trend.get("sources", [])
    assignees = trend.get("assignees", [])
    targets   = trend.get("targets", [])
    ws        = trend.get("whitespace", False)
    score     = trend.get("score", 0)

    vars_ = dict(
        cls           = cls,
        evidence      = evidence,
        src_str       = _src_str(sources),
        assignee_note = (f"{assignees[0]} leading" if assignees else "multiple groups"),
        vel_desc      = _vel_desc(score),
        whitespace_note  = _whitespace_note(ws, sources),
        assignee_detail  = _assignee_detail(assignees),
        targets_str   = _targets_str(targets),
    )

    t1 = _ONE_LINER_TEMPLATES.get(badge, _ONE_LINER_TEMPLATES["steady"])
    t2 = _DETAIL_TEMPLATES.get(badge, _DETAIL_TEMPLATES["steady"])
    t3 = _ACTION_TEMPLATES.get(badge, _ACTION_TEMPLATES["steady"])

    trend["one_liner"] = t1.safe_substitute(vars_)
    trend["detail"]    = t2.safe_substitute(vars_)
    trend["action"]    = t3.safe_substitute(vars_)
    return trend


# ── LLM-enhanced path (opt-in) ────────────────────────────────────────────────

_LLM_SYSTEM = """You are an editorial assistant for MolTrend, a medicinal chemistry intelligence tool.
Your job is to write punchy, precise card copy for building-block trend cards.
Rules:
- Use ONLY the numbers and facts provided in the JSON — never invent data.
- one_liner: 1 sentence, ≤15 words, present tense, concrete.
- detail: 2–3 sentences, specific, chemist-audience.
- action: 1 sentence, imperative, actionable.
- No filler phrases ("it is worth noting", "please be aware").
- Output ONLY valid JSON with keys: one_liner, detail, action.
"""


def write_copy_llm(trend: dict) -> dict:
    """
    Try to enhance copy with Claude Haiku. Falls back to deterministic on any error.
    Requires ANTHROPIC_API_KEY environment variable.
    """
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return write_copy_deterministic(trend)

    try:
        import anthropic
    except ImportError:
        logger.warning("anthropic package not installed — using deterministic copy")
        return write_copy_deterministic(trend)

    # First get deterministic copy as context
    trend = write_copy_deterministic(trend)

    prompt_data = {
        "bb_class":   trend.get("bb_class"),
        "badge":      trend.get("badge"),
        "score":      trend.get("score"),
        "evidence":   trend.get("evidence"),
        "sources":    trend.get("sources"),
        "assignees":  trend.get("assignees", [])[:3],
        "targets":    trend.get("targets", [])[:4],
        "whitespace": trend.get("whitespace"),
        "sparkline":  trend.get("sparkline"),
        "draft_one_liner": trend.get("one_liner"),
        "draft_detail":    trend.get("detail"),
        "draft_action":    trend.get("action"),
    }

    try:
        client = anthropic.Anthropic(api_key=api_key)
        msg = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=400,
            system=_LLM_SYSTEM,
            messages=[{
                "role": "user",
                "content": f"Rewrite copy for this trend card:\n{json.dumps(prompt_data, indent=2)}"
            }]
        )
        text = msg.content[0].text.strip()
        # Extract JSON
        start = text.find("{")
        end   = text.rfind("}") + 1
        if start >= 0 and end > start:
            out = json.loads(text[start:end])
            trend["one_liner"] = out.get("one_liner", trend["one_liner"])
            trend["detail"]    = out.get("detail",    trend["detail"])
            trend["action"]    = out.get("action",    trend["action"])
    except Exception as e:
        logger.warning("LLM editorial failed (%s) — keeping deterministic copy", e)

    return trend


def annotate_trends(trends: list[dict], use_llm: bool = True) -> list[dict]:
    """Write copy for all trend cards. use_llm=True attempts LLM if key available."""
    result = []
    for t in trends:
        if use_llm and os.environ.get("ANTHROPIC_API_KEY"):
            t = write_copy_llm(t)
        else:
            t = write_copy_deterministic(t)
        result.append(t)
    return result
