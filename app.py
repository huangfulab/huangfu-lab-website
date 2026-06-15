from flask import Flask, render_template, jsonify, request
import os
import re
import yaml
import glob
import time
import pandas as pd
from pathlib import Path
from perturbseq_bp import perturbseq_bp, _fmt_time_ago, _LAST_COMMIT_TS

app = Flask(__name__)
app.register_blueprint(perturbseq_bp)
DATA_DIR = Path(__file__).resolve().parent / "networks"
_cache = {}

# ── Lab homepage data (Google Sheets or YAML fallback) ──────────────────────

_LAB_YAML_SRC = Path("/data1/huangfud/torred1/sandbox/sandbox018-huangfu_lab_website/Huangfu-lab-website")
_CITATIONS_PATH = _LAB_YAML_SRC / "_data" / "citations.yaml"
_SHEETS_TTL = 300  # seconds
_sheets_cache = {}  # keys: "data", "fetched_at"

# Load the full publication list from ORCID-generated YAML once at startup.
with open(_CITATIONS_PATH) as _f:
    _ALL_CITATIONS = sorted(
        [c for c in (yaml.safe_load(_f) or []) if c.get("date")],
        key=lambda c: c["date"],
        reverse=True,
    )


def _safe_url(val):
    """Allow only http/https URLs; blank out javascript:, data:, etc."""
    v = str(val or "").strip()
    return v if v.startswith(("http://", "https://")) else ""


def _fetch_from_sheets():
    """Pull all editable content from Google Sheets."""
    import gspread
    from google.oauth2.service_account import Credentials

    creds = Credentials.from_service_account_file(
        os.environ["GOOGLE_SHEETS_CREDENTIALS"],
        scopes=["https://www.googleapis.com/auth/spreadsheets.readonly"],
    )
    gc = gspread.authorize(creds)
    sh = gc.open_by_key(os.environ["LAB_SHEET_ID"])

    # Description tab → {key: value} dict
    desc_rows = sh.worksheet("Description").get_all_records()
    description = {r["key"]: r["value"] for r in desc_rows if r.get("key")}

    # Members tab
    raw_members = sh.worksheet("Members").get_all_records()
    members = []
    for r in raw_members:
        if not r.get("name"):
            continue
        links = {}
        for field in ("email", "github", "orcid", "twitter"):
            if r.get(field):
                links[field] = r[field]
        if r.get("website"):
            links["website"] = _safe_url(r["website"])
        members.append({
            "name": r["name"],
            "role": r.get("role", ""),
            "affiliation": r.get("affiliation", ""),
            "bio": r.get("bio", ""),
            "image_url": _safe_url(r.get("image_url", "")),
            "links": links,
            "_order": int(r["order"]) if str(r.get("order", "")).isdigit() else 99,
        })
    members.sort(key=lambda x: (x["role"] != "principal-investigator", x["_order"], x["name"]))

    # News tab
    raw_news = sh.worksheet("News").get_all_records()
    news = sorted(
        [{**r, "link": _safe_url(r.get("link", ""))} for r in raw_news if r.get("title")],
        key=lambda r: r.get("date", ""),
        reverse=True,
    )

    # Featured Publications tab → set of DOI strings
    feat_rows = sh.worksheet("Featured Publications").get_all_records()
    featured_dois = {r["doi"].strip() for r in feat_rows if r.get("doi")}
    featured_citations = [c for c in _ALL_CITATIONS if c.get("id") in featured_dois]

    # Projects tab
    raw_projects = sh.worksheet("Projects").get_all_records()
    projects = []
    for r in raw_projects:
        if not r.get("title"):
            continue
        projects.append({
            "title": r["title"],
            "subtitle": r.get("subtitle", ""),
            "description": r.get("description", ""),
            "link": _safe_url(r.get("link", "")),
            "tags": [t.strip() for t in str(r.get("tags", "")).split(",") if t.strip()],
            "image_url": _safe_url(r.get("image_url", "")),
            "_order": int(r["order"]) if str(r.get("order", "")).isdigit() else 99,
        })
    projects.sort(key=lambda x: x["_order"])

    return {
        "description": description,
        "members": members,
        "news": news,
        "featured_citations": featured_citations,
        "projects": projects,
        "citations": _ALL_CITATIONS,
    }


def _load_yaml_fallback():
    """Fallback: read members/projects from Jekyll YAML (no Sheets credentials)."""
    raw_members = []
    for md_path in sorted(glob.glob(str(_LAB_YAML_SRC / "_members" / "*.md"))):
        with open(md_path) as f:
            content = f.read()
        m = re.match(r"^---\n(.*?)\n---\n?(.*)", content, re.DOTALL)
        if m:
            meta = yaml.safe_load(m.group(1)) or {}
            meta["bio"] = m.group(2).strip()
            meta.setdefault("image_url", "")
            meta.setdefault("_order", 99)
            raw_members.append(meta)
    raw_members.sort(key=lambda x: (x.get("role") != "principal-investigator", x.get("name", "")))

    with open(_LAB_YAML_SRC / "_data" / "projects.yaml") as f:
        raw_projects = yaml.safe_load(f) or []
    for p in raw_projects:
        p.setdefault("image_url", "")

    return {
        "description": {
            "description": (
                "The Huangfu Laboratory for Developmental and Stem Cell Biology is in "
                "the DevBio Program at Memorial Sloan Kettering Cancer Center, New York, NY."
            )
        },
        "members": raw_members,
        "news": [],
        "featured_citations": [],
        "projects": raw_projects,
        "citations": _ALL_CITATIONS,
    }


def _get_lab_data():
    """Return cached Sheets data, re-fetching if the TTL has expired."""
    if not (os.environ.get("GOOGLE_SHEETS_CREDENTIALS") and os.environ.get("LAB_SHEET_ID")):
        return _load_yaml_fallback()

    now = time.time()
    if _sheets_cache.get("data") and now - _sheets_cache.get("fetched_at", 0) < _SHEETS_TTL:
        return _sheets_cache["data"]

    try:
        data = _fetch_from_sheets()
        _sheets_cache["data"] = data
        _sheets_cache["fetched_at"] = now
        return data
    except Exception as e:
        app.logger.error("Sheets fetch failed: %s", e)
        return _sheets_cache.get("data") or _load_yaml_fallback()


# ── Network data ─────────────────────────────────────────────────────────────

def load_network(level):
    if level in _cache:
        return _cache[level]

    nodes_df = pd.read_csv(DATA_DIR / f"{level}_nodes.tsv", sep="\t")
    edges_df = pd.read_csv(DATA_DIR / f"{level}_edges.tsv", sep="\t")

    nodes = []
    for _, r in nodes_df.iterrows():
        node = {
            "id":            str(r["id"]),
            "n_genes":       int(r["n_genes"]),
            "color":         str(r["color"]),
            "within_mean_z": round(float(r["within_mean_z"]), 2),
        }
        if level == "submodule":
            node["supermodule"] = str(r["supermodule"])
        nodes.append(node)

    edges = []
    for _, r in edges_df.iterrows():
        mz = float(r["mean_z"])
        edges.append({
            "source":    str(r["source"]),
            "target":    str(r["target"]),
            "mean_z":    round(mz, 2),
            "abs_mean_z": round(abs(mz), 2),
        })

    _cache[level] = {"nodes": nodes, "edges": edges}
    return _cache[level]


# ── Routes ───────────────────────────────────────────────────────────────────

@app.route("/lab")
def lab():
    return render_template("lab/index.html")


@app.route("/lab_homepage_github")
def lab_homepage_github():
    data = _get_lab_data()
    return render_template("lab/github_homepage.html", **data)


@app.route("/api/lab/refresh", methods=["POST"])
def lab_refresh():
    token = os.environ.get("LAB_REFRESH_TOKEN")
    if token and request.headers.get("Authorization") != f"Bearer {token}":
        return jsonify({"error": "Unauthorized"}), 401
    _sheets_cache.clear()
    return jsonify({"ok": True})


@app.route("/")
def landing():
    return render_template("landing.html", last_updated=_fmt_time_ago(_LAST_COMMIT_TS))


@app.route("/modules")
def modules():
    return render_template("modules.html")


@app.route("/network")
def network():
    return render_template("index.html")


@app.route("/api/network/<level>")
def api_network(level):
    if level not in ("supermodule", "submodule"):
        return jsonify({"error": "invalid level"}), 400
    return jsonify(load_network(level))


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5050))
    app.run(debug=True, host="0.0.0.0", port=port)
