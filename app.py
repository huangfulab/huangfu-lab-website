from flask import Flask, render_template, jsonify
import os
import re
import yaml
import glob
import pandas as pd
from pathlib import Path
from perturbseq_bp import perturbseq_bp, _fmt_time_ago, _LAST_COMMIT_TS

app = Flask(__name__)
app.register_blueprint(perturbseq_bp)
DATA_DIR = Path(__file__).resolve().parent / "networks"
_cache = {}

LAB_JEKYLL_SRC = Path("/data1/huangfud/torred1/sandbox/sandbox018-huangfu_lab_website/Huangfu-lab-website")

def _load_lab_data():
    # Publications
    citations_path = LAB_JEKYLL_SRC / "_data" / "citations.yaml"
    with open(citations_path) as f:
        raw = yaml.safe_load(f) or []
    citations = sorted(
        [c for c in raw if c.get("date")],
        key=lambda c: c["date"],
        reverse=True,
    )

    # Team members
    members = []
    for md_path in sorted(glob.glob(str(LAB_JEKYLL_SRC / "_members" / "*.md"))):
        with open(md_path) as f:
            content = f.read()
        m = re.match(r"^---\n(.*?)\n---\n?(.*)", content, re.DOTALL)
        if m:
            meta = yaml.safe_load(m.group(1)) or {}
            meta["bio"] = m.group(2).strip()
            members.append(meta)
    # PI first, then rest
    members.sort(key=lambda x: (x.get("role") != "principal-investigator", x.get("name", "")))

    # Projects
    projects_path = LAB_JEKYLL_SRC / "_data" / "projects.yaml"
    with open(projects_path) as f:
        projects = yaml.safe_load(f) or []

    return citations, members, projects

CITATIONS, MEMBERS, PROJECTS = _load_lab_data()


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


@app.route("/lab")
def lab():
    return render_template("lab/index.html")


@app.route("/lab_homepage_github")
def lab_homepage_github():
    return render_template("lab/github_homepage.html",
                           citations=CITATIONS,
                           members=MEMBERS,
                           projects=PROJECTS)


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
