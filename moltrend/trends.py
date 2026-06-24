"""
trends.py — Scoring, sparklines, heat score, badges, PCA-2D embedding.

Key design:
  - Corpus-spanning buckets: split [min_year, max_year] into N=8 equal buckets.
    Recent live corpora → roughly yearly; historical → wider. Both give a
    meaningful sparkline.
  - Heat score is a transparent composite (weights documented below).
  - Whitespace = leading indicator: patent/paper signal with zero approvals.
  - PCA-via-SVD 2D embedding of Morgan fingerprints (no sklearn needed).
"""
from __future__ import annotations

import math
from collections import defaultdict
from datetime import datetime
from typing import Optional

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

# ── Constants ────────────────────────────────────────────────────────────────
N_BUCKETS = 8
MIN_SUPPORT = 3          # minimum molecules to surface a trend (raise for large corpora)
TOP_N_TRENDS = 12        # max trend cards surfaced

# Heat score weights (must sum to 1.0)
W_VELOCITY     = 0.34
W_ACCELERATION = 0.24
W_RECENT_VOL   = 0.24
W_NOVELTY      = 0.10
W_WHITESPACE   = 0.08

# Badge thresholds
HEAT_HOT    = 66
VEL_RISING  = 0.15

# Color palette (assigned by rank)
COLORS = [
    "#E8364E", "#0EA573", "#E67E22", "#6C5CE7",
    "#0891B2", "#DC6B30", "#7F8C9B", "#A29BFE",
    "#00B894", "#D63031", "#FDCB6E", "#2D3436",
]


# ── Date helpers ─────────────────────────────────────────────────────────────

def _parse_year(date_str: str) -> Optional[int]:
    """Extract year from 'YYYY-MM-DD' or 'YYYY'."""
    if not date_str:
        return None
    try:
        return int(str(date_str)[:4])
    except (ValueError, TypeError):
        return None


def _corpus_years(records: list[dict]) -> tuple[int, int]:
    """Return (min_year, max_year) from all records with a date."""
    years = [y for r in records if (y := _parse_year(r.get("date", ""))) is not None]
    if not years:
        now = datetime.now().year
        return now - 7, now
    return min(years), max(years)


def _bucket_index(year: int, min_year: int, max_year: int) -> int:
    """Map a year to a bucket index [0, N_BUCKETS-1]."""
    span = max(max_year - min_year, 1)
    raw = (year - min_year) / span * N_BUCKETS
    return min(int(raw), N_BUCKETS - 1)


def _make_sparkline(dates: list[str], min_year: int, max_year: int) -> list[int]:
    buckets = [0] * N_BUCKETS
    for d in dates:
        y = _parse_year(d)
        if y is not None:
            buckets[_bucket_index(y, min_year, max_year)] += 1
    return buckets


def _bucket_start_years(min_year: int, max_year: int) -> list[int]:
    span = max_year - min_year
    step = max(span / N_BUCKETS, 1)
    return [int(min_year + i * step) for i in range(N_BUCKETS)]


# ── Velocity & acceleration ───────────────────────────────────────────────────

def _velocity(spark: list[int]) -> float:
    """Slope of the last half vs first half, normalized by mean."""
    n = len(spark)
    if n < 2:
        return 0.0
    half = n // 2
    first = sum(spark[:half]) / max(half, 1)
    last  = sum(spark[half:]) / max(n - half, 1)
    mean  = (first + last) / 2 or 1
    return (last - first) / mean


def _acceleration(spark: list[int]) -> float:
    """Compare last third vs middle third."""
    n = len(spark)
    if n < 3:
        return 0.0
    t = n // 3
    mid  = sum(spark[t:2*t]) / max(t, 1)
    last = sum(spark[2*t:])  / max(n - 2*t, 1)
    mean = (mid + last) / 2 or 1
    return (last - mid) / mean


def _recent_volume(spark: list[int]) -> int:
    """Molecules in last two buckets."""
    return sum(spark[-2:])


# ── Normalization ─────────────────────────────────────────────────────────────

def _minmax(values: list[float]) -> list[float]:
    mn, mx = min(values), max(values)
    span = mx - mn or 1
    return [(v - mn) / span for v in values]


# ── Heat score & badge ────────────────────────────────────────────────────────

def _heat(vel_n: float, acc_n: float, rec_n: float, nov: float, ws: float) -> float:
    raw = (W_VELOCITY * vel_n + W_ACCELERATION * acc_n +
           W_RECENT_VOL * rec_n + W_NOVELTY * nov + W_WHITESPACE * ws)
    return round(raw * 100, 1)


def _badge(spark: list[int], heat: float, vel: float) -> str:
    recent_buckets = spark[-2:] if len(spark) >= 2 else spark
    older = spark[:-2] if len(spark) > 2 else []
    first_appeared_recently = all(v == 0 for v in older) and any(v > 0 for v in recent_buckets)
    if first_appeared_recently and vel > 0:
        return "emerging"
    if heat >= HEAT_HOT and vel > 0:
        return "hot"
    if vel > VEL_RISING:
        return "rising"
    return "steady"


def _trajectory(vel: float, acc: float) -> str:
    if vel > 0 and acc > 0:
        return "accelerating"
    if vel > 0:
        return "stable"
    return "emerging"


# ── PCA-2D embedding (Morgan fingerprints, SVD) ───────────────────────────────

def _morgan_fp_array(smiles_list: list[str], nbits: int = 512) -> Optional[list[list[int]]]:
    """Compute Morgan fingerprints for a list of SMILES. Returns None if RDKit absent."""
    if not RDKIT_OK:
        return None
    fps = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits)
                fps.append(list(fp))
            else:
                fps.append([0] * nbits)
        except Exception:
            fps.append([0] * nbits)
    return fps


def _pca2d_svd(matrix: list[list[float]]) -> list[tuple[float, float]]:
    """
    Pure-Python PCA via covariance SVD. Projects n×d matrix to n×2.
    Returns list of (x, y) in [0, 1] range.
    """
    n = len(matrix)
    if n == 0:
        return []
    d = len(matrix[0])

    # Center
    mean = [sum(matrix[i][j] for i in range(n)) / n for j in range(d)]
    X = [[matrix[i][j] - mean[j] for j in range(d)] for i in range(n)]

    # Gram matrix X @ X^T (n×n, cheaper when n << d)
    gram = [[sum(X[i][k] * X[jj][k] for k in range(d)) for jj in range(n)] for i in range(n)]

    # Power iteration for top-2 eigenvectors of gram
    def power_iter(M, iters=50):
        import random
        v = [random.gauss(0, 1) for _ in range(len(M))]
        for _ in range(iters):
            Mv = [sum(M[i][j] * v[j] for j in range(len(M))) for i in range(len(M))]
            norm = math.sqrt(sum(x*x for x in Mv)) or 1
            v = [x / norm for x in Mv]
        lam = sum(sum(M[i][j] * v[j] for j in range(len(M))) * v[i] for i in range(len(M)))
        return lam, v

    def deflate(M, lam, v):
        n = len(M)
        return [[M[i][j] - lam * v[i] * v[j] for j in range(n)] for i in range(n)]

    lam1, v1 = power_iter(gram)
    gram2 = deflate(gram, lam1, v1)
    _, v2 = power_iter(gram2)

    # Project X onto v1, v2 (these are in sample space; back-project via X^T)
    def project(v):
        proj = [sum(X[i][j] * v[i] for i in range(n)) for j in range(d)]
        norm = math.sqrt(sum(p*p for p in proj)) or 1
        return [p / norm for p in proj]

    pc1 = project(v1)
    pc2 = project(v2)

    coords = [(sum(X[i][j] * pc1[j] for j in range(d)),
               sum(X[i][j] * pc2[j] for j in range(d))) for i in range(n)]

    # Normalize to [0.05, 0.95]
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    xmn, xmx = min(xs), max(xs)
    ymn, ymx = min(ys), max(ys)
    xspan = (xmx - xmn) or 1
    yspan = (ymx - ymn) or 1

    return [(0.05 + 0.9 * (x - xmn) / xspan,
             0.05 + 0.9 * (y - ymn) / yspan) for x, y in coords]


# ── Main scoring function ─────────────────────────────────────────────────────

def score_groups(
    class_groups: dict[str, dict],
    brics_groups: dict[str, dict],
    records: list[dict],
    min_support: int = MIN_SUPPORT,
) -> tuple[list[dict], list[dict]]:
    """
    Score class groups and produce trend cards + compound points.

    Returns:
      (trends, compounds)
      trends:    list of trend dicts matching data.json schema
      compounds: list of compound dicts for t-SNE canvas
    """
    if not records:
        return [], []

    min_year, max_year = _corpus_years(records)
    bucket_years = _bucket_start_years(min_year, max_year)

    # ── Score each class group ────────────────────────────────────────────────
    raw_scores = []
    for cls_name, g in class_groups.items():
        evidence = len(g["members"])
        if evidence < min_support:
            continue

        spark = _make_sparkline(g["dates"], min_year, max_year)
        vel   = _velocity(spark)
        acc   = _acceleration(spark)
        rec   = _recent_volume(spark)
        nov   = 1.0 if max(spark[:-2], default=0) == 0 else 0.0
        ws    = 1.0 if (g["approval_count"] == 0 and (g["patent_count"] + g["paper_count"]) > 0) else 0.0

        raw_scores.append({
            "_cls":        cls_name,
            "_g":          g,
            "evidence":    evidence,
            "sparkline":   spark,
            "vel":         vel,
            "acc":         acc,
            "rec_vol":     rec,
            "novelty":     nov,
            "whitespace":  ws,
            "ws_flag":     ws > 0,
        })

    if not raw_scores:
        return [], []

    # Normalize each metric across all groups
    vels  = _minmax([r["vel"]    for r in raw_scores])
    accs  = _minmax([r["acc"]    for r in raw_scores])
    recs  = _minmax([r["rec_vol"] for r in raw_scores])

    # Compute heat scores
    for i, r in enumerate(raw_scores):
        r["heat"] = _heat(vels[i], accs[i], recs[i], r["novelty"], r["whitespace"])

    # Sort by heat, take top N
    raw_scores.sort(key=lambda x: x["heat"], reverse=True)
    raw_scores = raw_scores[:TOP_N_TRENDS]

    # ── Build trend cards ─────────────────────────────────────────────────────
    trends = []
    for rank, r in enumerate(raw_scores):
        cls_name = r["_cls"]
        g        = r["_g"]
        heat     = r["heat"]
        spark    = r["sparkline"]
        vel      = r["vel"]
        acc      = r["acc"]
        badge    = _badge(spark, heat, vel)
        traj     = _trajectory(vel, acc)
        color    = COLORS[rank % len(COLORS)]
        tid      = f"t{rank + 1}"

        # Source set
        src_list = sorted(g["sources"])

        # Assignees (top 5 by frequency)
        assignee_counts: dict[str, int] = defaultdict(int)
        for rec in records:
            if rec.get("assignee") and cls_name in classify_molecule_cached(rec.get("smiles", "")):
                assignee_counts[rec["assignee"]] += 1
        top_assignees = sorted(assignee_counts, key=lambda k: -assignee_counts[k])[:5]

        trends.append({
            "id":          tid,
            "headline":    _headline(cls_name, badge, traj),
            "bb_class":    cls_name,
            "score":       heat,
            "badge":       badge,
            "trajectory":  traj,
            "one_liner":   "",          # filled by editorial.py
            "detail":      "",
            "action":      "",
            "evidence":    r["evidence"],
            "sources":     src_list,
            "targets":     list(g["targets"])[:6],
            "assignees":   top_assignees,
            "whitespace":  r["ws_flag"],
            "suppliers":   {"enamine": 0, "molport": 0, "mcule": 0},  # enriched later
            "sparkline":   spark,
            "years":       bucket_years,
            "color":       color,
            "smiles":      g["member_smiles"][:3],
            "primary_smiles": g.get("primary_smiles", (g["member_smiles"] or [""])[0]),
        })

    # ── Build compound canvas points ──────────────────────────────────────────
    # Sample up to 200 compounds for the canvas (spread across trends)
    canvas_records = _sample_canvas_records(records, trends, max_per_trend=25)

    # Compute PCA-2D
    smiles_list = [r["smiles"] for r in canvas_records]
    fps = _morgan_fp_array(smiles_list)
    if fps and len(fps) == len(canvas_records):
        try:
            coords = _pca2d_svd(fps)
        except Exception:
            coords = [(0.5, 0.5)] * len(canvas_records)
    else:
        # Fallback: spread evenly
        coords = [(_pseudo_x(i, len(canvas_records)), _pseudo_y(i, len(canvas_records)))
                  for i in range(len(canvas_records))]

    # Map molecule → trend id
    trend_map = {smi: t["id"] for t in trends for smi in t["smiles"]}
    # Assign to nearest trend by structural class membership
    cls_to_tid = {t["bb_class"]: t["id"] for t in trends}

    compounds = []
    for i, rec in enumerate(canvas_records):
        smi = rec["smiles"]
        # Assign trend
        mol_classes = classify_molecule_cached(smi)
        tid = None
        for c in mol_classes:
            if c in cls_to_tid:
                tid = cls_to_tid[c]
                break
        if not tid:
            tid = trends[0]["id"] if trends else "t1"

        x, y = coords[i] if i < len(coords) else (0.5, 0.5)
        compounds.append({
            "smiles": smi,
            "name":   rec.get("name") or rec.get("ref_id") or smi[:20],
            "trend":  tid,
            "tsne":   [round(x, 4), round(y, 4)],
            "mw":     rec.get("mw", 0.0),
            "logd":   rec.get("clogp", 0.0),
            "hbd":    rec.get("hbd", 0),
            "hba":    rec.get("hba", 0),
            "tpsa":   rec.get("tpsa", 0.0),
            "rotb":   rec.get("rotb", 0),
        })

    return trends, compounds


# ── Helpers ───────────────────────────────────────────────────────────────────

# Cache for classify_molecule calls (to avoid re-running SMARTS on same SMILES)
_cls_cache: dict[str, list[str]] = {}


def classify_molecule_cached(smiles: str) -> list[str]:
    if smiles not in _cls_cache:
        from moltrend.classify import classify_molecule
        _cls_cache[smiles] = classify_molecule(smiles)
    return _cls_cache[smiles]


def _headline(cls_name: str, badge: str, traj: str) -> str:
    pfx = {"hot": "🔥", "rising": "📈", "emerging": "✨", "steady": "—"}.get(badge, "")
    return f"{cls_name} — {traj.capitalize()}"


def _sample_canvas_records(records: list[dict], trends: list[dict], max_per_trend: int) -> list[dict]:
    """Sample records so each trend class is represented."""
    cls_to_tid = {t["bb_class"]: t["id"] for t in trends}
    buckets: dict[str, list] = defaultdict(list)
    others: list = []

    for rec in records:
        classes = classify_molecule_cached(rec.get("smiles", ""))
        placed = False
        for c in classes:
            if c in cls_to_tid:
                buckets[c].append(rec)
                placed = True
                break
        if not placed:
            others.append(rec)

    sampled = []
    for cls, recs in buckets.items():
        sampled.extend(recs[:max_per_trend])
    sampled.extend(others[:max(0, 200 - len(sampled))])
    return sampled[:200]


def _pseudo_x(i: int, n: int) -> float:
    return 0.05 + 0.9 * (i % max(int(n ** 0.5), 1)) / max(int(n ** 0.5), 1)


def _pseudo_y(i: int, n: int) -> float:
    row = int(n ** 0.5)
    return 0.05 + 0.9 * (i // max(row, 1)) / max(row, 1)
