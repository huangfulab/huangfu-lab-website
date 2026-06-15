from flask import Flask, render_template, jsonify
import json
import os
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
def _clean_citation(c):
    if c.get("publisher"):
        c["publisher"] = html.unescape(c["publisher"])
    return c

with open(_CITATIONS_PATH) as _f:
    _ALL_CITATIONS = sorted(
        [_clean_citation(c) for c in (yaml.safe_load(_f) or []) if c.get("date")],
        key=lambda c: c["date"],
        reverse=True,
    )


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
    return render_template("lab/index.html", citations=_ALL_CITATIONS[:3])


@app.route("/research")
def lab_research():
    return render_template("lab/research.html", current_page="research")


@app.route("/publications")
def lab_publications():
    return render_template("lab/publications.html", current_page="publications", citations=_ALL_CITATIONS)


@app.route("/announcements")
def lab_announcements():
    return render_template("lab/announcements.html", current_page="announcements")


_TEAM_PATH = Path(__file__).resolve().parent / "static" / "team.json"
with open(_TEAM_PATH) as _tf:
    _TEAM = json.load(_tf)  # flat list; each entry has alumni: true/false


@app.route("/team")
def lab_team():
    return render_template("lab/team.html", current_page="team", team=_TEAM)


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


@app.route("/sitemap.xml")
def sitemap_xml():
    urls = [
        "/", "/research", "/publications", "/announcements",
        "/team", "/resources", "/contact",
        "/tf-perturbseq", "/modules", "/network",
    ]
    xml = render_template("sitemap.xml", urls=urls)
    return xml, 200, {"Content-Type": "application/xml"}


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5050))
    app.run(debug=True, host="0.0.0.0", port=port)
