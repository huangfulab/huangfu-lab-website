from flask import Flask, render_template, jsonify
import json
import os
import sqlite3
import yaml
import html
import pandas as pd
from pathlib import Path
from perturbseq_bp import perturbseq_bp

app = Flask(__name__)
app.register_blueprint(perturbseq_bp)
DATA_DIR = Path(__file__).resolve().parent / "networks"
_cache = {}

_CITATIONS_PATH = Path(__file__).resolve().parent / "citations.yaml"
_PREPRINT_PUBLISHERS = {
    "10.1101/": "bioRxiv",
    "10.64898/": "openRxiv",
}

def _clean_citation(c):
    if c.get("publisher"):
        c["publisher"] = html.unescape(c["publisher"])
    doi = c.get("id", "").removeprefix("doi:")
    for prefix, name in _PREPRINT_PUBLISHERS.items():
        if doi.startswith(prefix):
            c["publisher"] = name
            break
    return c

with open(_CITATIONS_PATH) as _f:
    _ALL_CITATIONS = sorted(
        [_clean_citation(c) for c in (yaml.safe_load(_f) or []) if c.get("date")],
        key=lambda c: c["date"],
        reverse=True,
    )

_PREPRINT_NAMES = set(_PREPRINT_PUBLISHERS.values())

_LAST_AUTHOR_CITATIONS = [
    c for c in _ALL_CITATIONS
    if c.get("authors") and c["authors"][-1] == "Danwei Huangfu"
    and c.get("publisher") not in _PREPRINT_NAMES
]


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


@app.route("/")
@app.route("/lab")
def lab():
    return render_template("lab/index.html", citations=_LAST_AUTHOR_CITATIONS[:5], announcements=_load_announcements()[:5])


@app.route("/research")
def lab_research():
    return render_template("lab/research.html", current_page="research")


@app.route("/publications")
def lab_publications():
    return render_template("lab/publications.html", current_page="publications", citations=_ALL_CITATIONS)


@app.route("/announcements")
def lab_announcements():
    return render_template("lab/announcements.html", current_page="announcements", announcements=_load_announcements())


_TEAM_PATH = Path(__file__).resolve().parent / "static" / "json" / "team.json"
_ANNOUNCEMENTS_PATH = Path(__file__).resolve().parent / "static" / "json" / "announcements.json"

def _load_team():
    with open(_TEAM_PATH) as f:
        return json.load(f)

def _load_announcements():
    with open(_ANNOUNCEMENTS_PATH) as f:
        return json.load(f)


@app.route("/team")
def lab_team():
    return render_template("lab/team.html", current_page="team", team=_load_team())


@app.route("/tf-perturbseq")
def tf_perturbseq():
    return render_template("landing.html")


@app.route("/resources")
def lab_resources():
    return render_template("lab/resources.html", current_page="resources")


@app.route("/contact")
def lab_contact():
    return render_template("lab/contact.html", current_page="contact")


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


@app.route("/robots.txt")
def robots_txt():
    from flask import request as _req
    sitemap_url = _req.host_url.rstrip("/") + "/sitemap.xml"
    body = f"User-agent: *\nAllow: /\nSitemap: {sitemap_url}\n"
    return body, 200, {"Content-Type": "text/plain"}


_SITEMAP_DB_PATH = Path(__file__).resolve().parent / "data" / "db" / "tf-perturbseq-v5.db"

@app.route("/sitemap.xml")
def sitemap_xml():
    urls = [
        "/", "/research", "/publications", "/announcements",
        "/team", "/resources", "/contact",
        "/tf-perturbseq", "/tf-perturbseq/all-modules", "/modules", "/network",
    ]
    try:
        conn = sqlite3.connect(str(_SITEMAP_DB_PATH))
        genes = [r[0] for r in conn.execute(
            "SELECT gene_name FROM gene_table ORDER BY gene_name"
        ).fetchall()]
        urls.extend(f"/tf-perturbseq/gene/{g}" for g in genes)
        supermodules = [r[0] for r in conn.execute(
            "SELECT module_name FROM module_table "
            "WHERE source='hotspot_supermodule' AND module_name != 'unassigned' ORDER BY module_name"
        ).fetchall()]
        urls.extend(f"/tf-perturbseq/module/{m}" for m in supermodules)
        submodules = [r[0] for r in conn.execute(
            "SELECT module_name FROM module_table "
            "WHERE source='hotspot_submodule' AND module_name != 'unassigned' ORDER BY module_name"
        ).fetchall()]
        urls.extend(f"/tf-perturbseq/module/{m}" for m in submodules)
        gene_clusters = [r[0] for r in conn.execute(
            "SELECT module_name FROM module_table "
            "WHERE source='mfuzz_k7' AND module_name != 'cluster_7' ORDER BY module_name"
        ).fetchall()]
        urls.extend(f"/tf-perturbseq/module/GC{m.split('_')[1]}" for m in gene_clusters)
        conn.close()
    except Exception:
        pass
    xml = render_template("sitemap.xml", urls=urls)
    return xml, 200, {"Content-Type": "application/xml"}


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5050))
    app.run(debug=True, host="0.0.0.0", port=port)
