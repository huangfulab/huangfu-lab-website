#!/usr/bin/env python
"""Test suite for the public API at /endoderm-perturbseq/api/v1.

    python scripts/test_api_v1.py                 # in-process, against the real DB
    python scripts/test_api_v1.py --base-url URL  # smoke-test a running server
    python scripts/test_api_v1.py --group objects
    python scripts/test_api_v1.py --only aliasing --verbose

Deliberately dependency-free: the repo has no test tooling, and adding pytest to
requirements.txt would install it into the production venv on the next deploy.

Most cases are generated from blueprints/api_v1/catalog.py, so a new endpoint
gets covered automatically and an undocumented one fails the coverage case.
"""
import argparse
import ast
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(ROOT / "scripts"))

from blueprints.api_v1.catalog import ENDPOINTS, ENUMS  # noqa: E402
from blueprints.api_v1.core import API_V1_PREFIX  # noqa: E402
from blueprints.api_v1.examples import EXAMPLE_RESPONSES  # noqa: E402
from blueprints.perturbseq_bp import DB_PATH  # noqa: E402

SENTINEL = "__NO_SUCH_THING__"
CASES = []


def case(name, group="misc"):
    def deco(fn):
        CASES.append((group, name, fn))
        return fn
    return deco


def expect(cond, msg):
    if not cond:
        raise AssertionError(msg)


def _no_json_constants(c):
    raise AssertionError(f"response contains the non-JSON constant {c!r}")


class Resp:
    __slots__ = ("status", "mimetype", "text", "headers")

    def __init__(self, status, mimetype, text, headers):
        self.status, self.mimetype, self.text, self.headers = status, mimetype, text, headers

    @property
    def json(self):
        # parse_constant fires on NaN / Infinity / -Infinity, which json.dumps
        # would only emit if something bypassed the serializer.
        return json.loads(self.text, parse_constant=_no_json_constants)


class LocalClient:
    """Flask test client. url_map and the app object are available, so the
    structural cases can run."""

    def __init__(self):
        import app as appmod
        self.app = appmod.app
        self.client = self.app.test_client()

    def request(self, url, method="GET"):
        r = self.client.open(url, method=method)
        return Resp(r.status_code, r.mimetype, r.get_data(as_text=True), dict(r.headers))


class RemoteClient:
    """urllib against a live server, for post-deploy smoke tests."""

    structural = False

    def __init__(self, base):
        self.base = base.rstrip("/")

    def request(self, url, method="GET"):
        req = urllib.request.Request(self.base + url, method=method)
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                body = r.read().decode()
                return Resp(r.status, r.headers.get_content_type(), body, dict(r.headers))
        except urllib.error.HTTPError as e:
            body = e.read().decode()
            return Resp(e.code, e.headers.get_content_type(), body, dict(e.headers))


def get(client, url):
    return client.request(url)


# ── helpers over the catalog ─────────────────────────────────────────────────

def _relative(url):
    return url[len(API_V1_PREFIX):] if url.startswith(API_V1_PREFIX) else url


def with_query(url, fragment):
    """Set a query parameter, replacing any existing value for that key.

    Several example_urls carry parameters of their own. Appending blindly would
    produce '?q=SOX?page=abc' (malformed) or '?per_page=5&per_page=abc', where
    Werkzeug returns the first value and the malformed one is never tested.
    """
    base, _, existing = url.partition("?")
    key = fragment.split("=", 1)[0]
    kept = [kv for kv in existing.split("&") if kv and kv.split("=", 1)[0] != key]
    return base + "?" + "&".join(kept + [fragment])


def not_found_url(ep):
    """example_url with every path param replaced by something that cannot exist.

    Type-aware, so a numeric segment stays numeric — otherwise the rule does not
    match at all and the case tests the catch-all instead of the endpoint.
    """
    if not ep["path_params"]:
        return None
    path = ep["path"]
    for p in ep["path_params"]:
        if p["type"] == "integer":
            filler = "999999999"
        elif p.get("enum"):
            filler = ENUMS[p["enum"]][0]
        else:
            filler = SENTINEL
        path = path.replace("{" + p["name"] + "}", filler)
    return API_V1_PREFIX + path


def int_params(ep):
    return [p for p in ep["query_params"] if p["type"] == "integer"]


def float_params(ep):
    return [p for p in ep["query_params"] if p["type"] == "number"]


def enum_params(ep):
    return [p for p in ep["query_params"] if p.get("enum")]


def assert_envelope(body, ep):
    expect("data" in body, f"no 'data' key; got {sorted(body)}")
    if ep["returns"] == "entity":
        expect(isinstance(body["data"], dict), "entity 'data' must be an object")
        expect("links" in body, "entity response must carry 'links'")
        expect(body["links"].get("self"), "links.self must be set")
    elif ep["returns"] == "collection":
        expect(isinstance(body["data"], list), "collection 'data' must be a list")
        for key in ("page", "total", "per_page", "links"):
            expect(key in body, f"collection response missing {key!r}")


def assert_json_error(resp, url, expected_status=None):
    expect(resp.mimetype == "application/json",
           f"{url} returned {resp.mimetype}, expected application/json")
    if expected_status:
        expect(resp.status == expected_status,
               f"{url} returned {resp.status}, expected {expected_status}")
    body = resp.json
    for key in ("error", "status", "code"):
        expect(key in body, f"{url} error body missing {key!r}: {body}")
    blob = resp.text
    for leak in ("Traceback", "SELECT ", "sqlite3.", DB_PATH):
        expect(leak not in blob, f"{url} error body leaks {leak!r}")
    return body


# ── generated cases ──────────────────────────────────────────────────────────

def register_generated():
    for ep in ENDPOINTS:
        eid, group, url = ep["id"], ep["group"], ep["example_url"]

        def happy(client, ep=ep, url=url):
            r = get(client, url)
            expect(r.status == 200, f"{url} returned {r.status}")
            expect(r.mimetype == "application/json", f"{url} returned {r.mimetype}")
            body = r.json
            if ep["returns"] in ("entity", "collection"):
                assert_envelope(body, ep)
        case(f"{eid}: happy path", group)(happy)

        nf = not_found_url(ep)
        if nf:
            def missing(client, nf=nf):
                assert_json_error(get(client, nf), nf, 404)
            case(f"{eid}: unknown identifier is a JSON 404", group)(missing)

        for p in int_params(ep) + float_params(ep):
            def bad_number(client, url=url, p=p):
                for bad in (f"{p['name']}=abc", f"{p['name']}="):
                    u = with_query(url, bad)
                    body = assert_json_error(get(client, u), u, 400)
                    expect(body["code"] == "invalid_param",
                           f"{u} gave code {body['code']!r}, expected invalid_param")
                if "min" in p:
                    u = with_query(url, f"{p['name']}={p['min'] - 1}")
                    assert_json_error(get(client, u), u, 400)
                if "max" in p:
                    u = with_query(url, f"{p['name']}={p['max'] + 1}")
                    body = assert_json_error(get(client, u), u, 400)
                    expect(str(p["max"]) in body["error"],
                           f"{u} error should name the cap {p['max']}: {body['error']}")
            case(f"{eid}: {p['name']} rejects malformed and out-of-range", group)(bad_number)

        for p in enum_params(ep):
            if p in ep["path_params"]:
                continue
            def bad_enum(client, url=url, p=p):
                u = with_query(url, f"{p['name']}=__nope__")
                body = assert_json_error(get(client, u), u, 400)
                allowed = ENUMS[p["enum"]]
                expect(any(v in body["error"] for v in allowed),
                       f"{u} error should list the valid values: {body['error']}")
            case(f"{eid}: {p['name']} rejects an unknown value", group)(bad_enum)

        def unknown_arg(client, url=url):
            u = with_query(url, "zzz_not_a_param=1")
            body = assert_json_error(get(client, u), u, 400)
            expect(body["code"] == "unknown_param",
                   f"{u} gave code {body['code']!r}, expected unknown_param")
        case(f"{eid}: unknown query parameter is rejected", group)(unknown_arg)

        def wrong_method(client, url=url):
            r = client.request(url, method="POST")
            assert_json_error(r, url, 405)
        case(f"{eid}: POST returns a JSON 405", group)(wrong_method)

        def cors(client, url=url):
            r = get(client, url)
            expect(r.headers.get("Access-Control-Allow-Origin") == "*",
                   f"{url} missing CORS header, got {r.headers.get('Access-Control-Allow-Origin')!r}")
        case(f"{eid}: sends an open CORS header", group)(cors)


# ── hand-written cases: the aliasing rules ───────────────────────────────────

@case("all four gene-cluster aliases resolve to one module", "aliasing")
def _aliases(client):
    ids = {}
    for name in ("GC1", "gc1", "TC-1", "gene_cluster_1", "cluster_1"):
        r = get(client, f"{API_V1_PREFIX}/module/{name}")
        expect(r.status == 200, f"/module/{name} returned {r.status}")
        d = r.json["data"]
        ids[name] = d["module_id"]
        expect(d["id"] == "GC1", f"/module/{name} reported id {d['id']!r}, expected 'GC1'")
        expect(d["module_name"] == "cluster_1", f"/module/{name} lost the DB name")
    expect(len(set(ids.values())) == 1, f"aliases resolved to different modules: {ids}")


@case("ambiguous module name requires ?source=", "aliasing")
def _ambiguous(client):
    u = f"{API_V1_PREFIX}/module/unassigned"
    body = assert_json_error(get(client, u), u, 400)
    expect(body["code"] == "ambiguous_name", f"got code {body['code']!r}")
    expect("hotspot_submodule" in body["error"] and "hotspot_supermodule" in body["error"],
           f"error should name both collections: {body['error']}")

    seen = {}
    for source in ("hotspot_supermodule", "hotspot_submodule"):
        r = get(client, f"{u}?source={source}")
        expect(r.status == 200, f"?source={source} returned {r.status}")
        seen[source] = r.json["data"]["module_id"]
    expect(seen["hotspot_supermodule"] != seen["hotspot_submodule"],
           f"both sources resolved to the same module: {seen}")


@case("hidden gene cluster GC7 is not served", "aliasing")
def _hidden(client):
    u = f"{API_V1_PREFIX}/module/GC7"
    assert_json_error(get(client, u), u, 404)


@case("supermodule and submodule are told apart without ?source=", "aliasing")
def _super_vs_sub(client):
    sup = get(client, f"{API_V1_PREFIX}/module/DE-1").json["data"]
    sub = get(client, f"{API_V1_PREFIX}/module/DE-1.2").json["data"]
    expect(sup["source"] == "hotspot_supermodule", f"DE-1 resolved to {sup['source']}")
    expect(sub["source"] == "hotspot_submodule", f"DE-1.2 resolved to {sub['source']}")
    expect(sub["parent_module"] == "DE-1", f"DE-1.2 parent was {sub['parent_module']!r}")
    expect("DE-1.2" in sup["child_submodules"], "DE-1 does not list DE-1.2 as a child")


@case("mfuzz perturbation GSEA uses the right gene_set_collection", "aliasing")
def _mfuzz_collection(client):
    # gsea_tf_table holds 21 near-identical mfuzz collections; only
    # ..._top1000_DE is the canonical one. Picking another, or failing to map
    # mfuzz_k7 at all, yields a silently empty perturbation arm.
    # GC2 is used deliberately: cluster_1 genuinely has no perturbation rows,
    # so GC1 cannot distinguish "correct" from "broken".
    d = get(client, f"{API_V1_PREFIX}/module/GC2?include=tf_regulators").json["data"]
    expect(d["tf_regulators"], "GC2 returned no TF regulators at all")
    perturbed = {r["tf_gene_name"] for r in d["tf_regulators"]
                 if r["evidence"] in ("perturbation", "both")}
    # These are GC2's regulators under the canonical top1000_DE collection. The
    # other 20 mfuzz collections give a different, smaller set, so merely
    # asserting "some perturbation evidence" would not detect the wrong one.
    expected = {"ARID1A", "FOXH1", "SMARCC1", "SOX17"}
    expect(expected <= perturbed,
           f"GC2 perturbation regulators are {sorted(perturbed)}, expected to include "
           f"{sorted(expected)} — names.GSEA_COLLECTION is probably pointing at the wrong "
           f"mfuzz gene_set_collection (there are 21 near-identical ones)")


@case("mfuzz binding enrichment crosses the mfuzz/mfuzz_k7 boundary", "aliasing")
def _mfuzz_binding(client):
    # tf_module_enrichment labels this collection 'mfuzz' while module_table
    # calls it 'mfuzz_k7'.
    d = get(client, f"{API_V1_PREFIX}/module/GC1?include=tf_regulators").json["data"]
    expect(any(r["evidence"] in ("binding", "both") for r in d["tf_regulators"]),
           "GC1 has no binding-supported regulators")


@case("a gene synonym resolves and is reported in meta", "aliasing")
def _synonym(client):
    r = get(client, f"{API_V1_PREFIX}/gene/HNF3B")
    expect(r.status == 200, f"synonym lookup returned {r.status}")
    body = r.json
    expect(body["data"]["gene_name"] == "FOXA2", f"HNF3B resolved to {body['data']['gene_name']}")
    meta = body.get("meta", {})
    expect(meta.get("resolved_from") == "HNF3B", f"meta.resolved_from missing: {meta}")
    expect(meta.get("match_type") == "synonym", f"meta.match_type was {meta.get('match_type')!r}")


@case("a duplicated symbol resolves to its primary row", "aliasing")
def _primary_gene(client):
    # 129 gene_table rows share a symbol with a primary entry. Without the
    # primary_gene=1 preference the pick is whichever row SQLite reaches first,
    # so the same request can return different genes across DB rebuilds.
    r = get(client, f"{API_V1_PREFIX}/gene/Y_RNA")
    expect(r.status == 200, f"/gene/Y_RNA returned {r.status}")
    d = r.json["data"]
    expect(d["is_primary_symbol"] is True,
           f"Y_RNA resolved to a non-primary row (gene_id {d['gene_id']}, "
           f"{d['ensg_id']}) — the primary_gene=1 preference is not being applied")


@case("a canonical symbol reports no resolution metadata", "aliasing")
def _canonical(client):
    body = get(client, f"{API_V1_PREFIX}/gene/FOXA2").json
    expect("resolved_from" not in body.get("meta", {}),
           "canonical lookup should not claim it resolved from something else")


# ── hand-written cases: data contract ────────────────────────────────────────

@case("infinite odds ratios are null with a companion flag", "contract")
def _infinite_odds(client):
    d = get(client, f"{API_V1_PREFIX}/tf/ARID1A?level=all").json["data"]
    flagged = [r for r in d["module_regulation"] if r["odds_ratio_infinite"]]
    expect(flagged, "no infinite odds ratios found — expected some in tf_module_enrichment")
    for r in flagged:
        expect(r["odds_ratio"] is None,
               f"{r['module']} flagged infinite but odds_ratio is {r['odds_ratio']!r}")


@case("a TF that is neither perturbed nor bound is rejected as a TF", "contract")
def _not_a_tf(client):
    genes = get(client, f"{API_V1_PREFIX}/gene/Y_RNA")
    if genes.status != 200:
        return  # symbol absent in this build; nothing to assert
    r = get(client, f"{API_V1_PREFIX}/tf/Y_RNA")
    if r.status == 200:
        d = r.json["data"]
        expect(d["is_perturbed_tf"] or d["is_binding_tf"],
               "a 200 from /tf must mean at least one kind of TF evidence")
    else:
        body = assert_json_error(r, "/tf/Y_RNA", 404)
        expect(body["code"] in ("tf_not_found", "gene_not_found"), f"code was {body['code']!r}")


@case("module evidence merging tags perturbation, binding and both", "contract")
def _evidence(client):
    d = get(client, f"{API_V1_PREFIX}/tf/ARID1A").json["data"]
    kinds = {r["evidence"] for r in d["module_regulation"]}
    expect(kinds <= {"perturbation", "binding", "both"}, f"unexpected evidence values: {kinds}")
    expect("both" in kinds, "expected at least one module with both kinds of evidence")


@case("gene clusters never leak the internal cluster_N name as the id", "contract")
def _no_internal_names(client):
    d = get(client, f"{API_V1_PREFIX}/gene/SOX17").json["data"]
    for m in d["modules"]:
        if m["source"] == "mfuzz_k7":
            expect(m["display_name"].startswith("GC"),
                   f"mfuzz module surfaced as {m['display_name']!r}")


@case("bounded lists carry a truncation flag", "contract")
def _truncation(client):
    # (path, {list field: its documented cap})
    checks = [
        ("/gene/SOX17?include=go_terms,perturbation_effects,elements",
         {"go_terms": 500, "perturbation_effects": 500, "elements": 200}),
        ("/tf/ARID1A?level=all", {"module_regulation": 500}),
        ("/module/DE-1?include=enrichment,tf_regulators", {"enrichment": 200, "tf_regulators": 500}),
    ]
    for path, caps in checks:
        d = get(client, API_V1_PREFIX + path).json["data"]
        for key, cap in caps.items():
            flag = f"{key}_truncated"
            expect(flag in d, f"{path}: {key} has no companion {flag}")
            expect(isinstance(d[flag], bool), f"{path}: {flag} is not a boolean")
            expect(len(d[key]) <= cap,
                   f"{path}: {key} returned {len(d[key])} rows, over its {cap} cap")
            # The flag must track reality in both directions, or it is decoration.
            expect(d[flag] == (len(d[key]) >= cap),
                   f"{path}: {flag} is {d[flag]} but {key} has {len(d[key])} of {cap} rows")


@case("expensive blocks are opt-in", "contract")
def _opt_in(client):
    plain = get(client, f"{API_V1_PREFIX}/gene/SOX17").json["data"]
    expect("coexpression" not in plain, "coexpression should not be in the default payload")
    withit = get(client, f"{API_V1_PREFIX}/gene/SOX17?include=coexpression").json["data"]
    expect(withit.get("coexpression"), "include=coexpression returned nothing")
    expect(len(withit["coexpression"]) <= 500, "coexpression exceeded its cap")


@case("link/gene-module refuses an unfiltered dump", "contract")
def _filter_required(client):
    u = f"{API_V1_PREFIX}/link/gene-module"
    body = assert_json_error(get(client, u), u, 400)
    expect(body["code"] == "filter_required", f"code was {body['code']!r}")
    ok = get(client, f"{u}?module=DE-1&per_page=3")
    expect(ok.status == 200, f"filtered request returned {ok.status}")
    expect(ok.json["total"] > 0, "DE-1 returned no members")


@case("tf-gene binding is found through the composite key", "contract")
def _tf_gene_binding(client):
    # This is the CROSS JOIN case: without it the planner probes
    # idx_tf_bs_gene_id and reads a 184M-row table instead of the PK.
    d = get(client, f"{API_V1_PREFIX}/link/tf-gene/FOXA2/SOX17?include=datasets").json["data"]
    expect(d["binding_evidence"] is True, "FOXA2->SOX17 should have binding evidence")
    expect(len(d["datasets"]) == d["n_datasets_with_binding"],
           "n_datasets_with_binding disagrees with the dataset list")
    expect(d["n_datasets_with_binding"] > 0, "no binding datasets returned")
    expect(d["max_binding_score_A"] is not None, "no binding score returned")
    for ds in d["datasets"]:
        expect("binding_score_A" in ds, "dataset row is missing its binding score")


@case("tf-gene element lookup is gated on binding evidence", "contract")
def _tf_gene_gate(client):
    # The expensive peak-level query must not run when no dataset places the TF
    # at the gene — that gate is what keeps this endpoint affordable.
    d = get(client, f"{API_V1_PREFIX}/link/tf-gene/SOX17/A2ML1?include=elements").json["data"]
    expect(d["binding_evidence"] is False, "expected no binding for SOX17->A2ML1")
    expect(d["elements"] == [], f"elements should be empty, got {len(d['elements'])}")
    expect(d.get("elements_note"), "a skipped element lookup should say so")

    with_binding = get(
        client, f"{API_V1_PREFIX}/link/tf-gene/SOX17/ACTB?include=elements").json["data"]
    expect(with_binding["binding_evidence"] is True, "expected binding for SOX17->ACTB")
    expect(with_binding["elements"], "a TF-gene pair with binding should return elements")


@case("search escapes LIKE wildcards instead of interpreting them", "contract")
def _search_wildcards(client):
    # Unescaped, '%' matches every row and turns the query into a full scan.
    for q in ("%%", "__", "%_%"):
        r = get(client, f"{API_V1_PREFIX}/search?q={q}")
        expect(r.status == 200, f"search q={q!r} returned {r.status}")
        expect(r.json["data"] == [],
               f"search q={q!r} returned {len(r.json['data'])} results — "
               f"wildcards are being interpreted, not escaped")


@case("search prefix-matches names and reports what matched", "contract")
def _search_prefix(client):
    rows = get(client, f"{API_V1_PREFIX}/search?q=SOX").json["data"]
    expect(rows, "search for SOX returned nothing")
    genes = [r for r in rows if r["type"] in ("gene", "tf")]
    expect(genes, "search for SOX returned no genes")
    for r in genes:
        expect(r["name"].upper().startswith("SOX"),
               f"{r['name']!r} is not a SOX prefix match")
        expect(r["matched_on"], "every result must report matched_on")

    syn = get(client, f"{API_V1_PREFIX}/search?q=HNF3B").json["data"]
    expect(any(r["matched_on"] == "synonym" for r in syn),
           "HNF3B should match as a synonym")


@case("tfs returns different populations per set", "contract")
def _tf_sets(client):
    totals = {}
    for name in ("perturbed", "binding", "all"):
        body = get(client, f"{API_V1_PREFIX}/tfs?set={name}&per_page=1").json
        totals[name] = body["total"]
    expect(totals["perturbed"] < totals["binding"],
           f"expected fewer perturbed than binding TFs, got {totals}")
    expect(totals["all"] >= totals["binding"],
           f"'all' must be the union, got {totals}")
    rows = get(client, f"{API_V1_PREFIX}/tfs?set=perturbed&per_page=500").json["data"]
    expect(all(r["is_perturbed_tf"] for r in rows), "set=perturbed returned a non-perturbed TF")


@case("an ATAC peak resolves identically by id and by name", "contract")
def _peak_by_name(client):
    by_id = get(client, f"{API_V1_PREFIX}/atac-peak/58938").json["data"]
    by_name = get(client, f"{API_V1_PREFIX}/atac-peak/{by_id['atac_peak_name']}").json["data"]
    expect(by_id["atac_peak_id"] == by_name["atac_peak_id"], "id and name disagree")
    expect(by_id["n_tf_peak_overlaps"] > 0, "expected TF overlaps on this peak")
    expect("accessibility" in by_id and by_id["accessibility"], "no accessibility profile")
    for tp in by_id["accessibility"]:
        expect("timepoint" in tp and "mean" in tp, f"malformed accessibility row: {tp}")


@case("the peak TF list is a separate, paginated endpoint", "contract")
def _peak_tfs(client):
    body = get(client, f"{API_V1_PREFIX}/atac-peak/2/tfs?per_page=5").json
    expect(len(body["data"]) <= 5, "per_page not honoured")
    expect(body["total"] is None,
           "total should be null here — counting distinct TFs is as costly as the query")
    for r in body["data"]:
        expect("tf_gene_name" in r and "n_datasets" in r, f"malformed row: {r}")


@case("GO accession forms are interchangeable", "contract")
def _go_forms(client):
    ids = set()
    for form in ("GO:0030183", "GO_0030183", "0030183"):
        r = get(client, f"{API_V1_PREFIX}/go-term/{form}")
        expect(r.status == 200, f"/go-term/{form} returned {r.status}")
        d = r.json["data"]
        ids.add(d["go_id"])
        expect(d["id"] == "GO:0030183", f"{form} reported id {d['id']!r}")
    expect(len(ids) == 1, f"accession forms resolved to different terms: {ids}")


@case("the two dataset namespaces are independent", "contract")
def _dataset_namespaces(client):
    tf = get(client, f"{API_V1_PREFIX}/dataset/tf/1").json["data"]
    ptm = get(client, f"{API_V1_PREFIX}/dataset/ptm/1").json["data"]
    expect(tf["namespace"] == "tf" and ptm["namespace"] == "ptm", "namespace not reported")
    expect(tf["dataset"] != ptm["dataset"],
           "id 1 returned the same record in both namespaces — they are being conflated")
    expect("n_tf_peaks" in tf, "TF dataset should report its peak count")
    expect("n_atac_peaks" in ptm, "PTM dataset should report its peak count")


@case("edges defaults to supermodules and honours the evidence filter", "contract")
def _edges(client):
    body = get(client, f"{API_V1_PREFIX}/edges?per_page=5").json
    expect(body["meta"]["level"] == "hotspot_supermodule", "unexpected default level")
    expect(body["total"] > 0, "no edges at the default level")
    for r in body["data"]:
        expect(r["evidence"] in ("perturbation", "binding", "both"), f"bad evidence: {r}")

    both = get(client, f"{API_V1_PREFIX}/edges?evidence=both&per_page=100").json
    expect(all(r["evidence"] == "both" for r in both["data"]),
           "evidence=both returned single-evidence edges")
    expect(both["total"] < body["total"], "evidence=both should be a strict subset")


@case("chromosome filters accept both 'chr7' and '7'", "contract")
def _chr_norm(client):
    a = get(client, f"{API_V1_PREFIX}/genes?chr=chr7&per_page=1").json["total"]
    b = get(client, f"{API_V1_PREFIX}/genes?chr=7&per_page=1").json["total"]
    expect(a == b and a > 0, f"chr7 gave {a}, 7 gave {b}")


@case("paging does not repeat or drop rows", "contract")
def _paging(client):
    p1 = get(client, f"{API_V1_PREFIX}/genes?per_page=25&page=1").json
    p2 = get(client, f"{API_V1_PREFIX}/genes?per_page=25&page=2").json
    ids1 = [r["gene_id"] for r in p1["data"]]
    ids2 = [r["gene_id"] for r in p2["data"]]
    expect(len(ids1) == 25 and len(ids2) == 25, "short page")
    expect(not set(ids1) & set(ids2), "page 1 and page 2 overlap")
    expect(p1["links"]["next"], "page 1 should advertise a next link")
    expect(p1["links"]["prev"] is None, "page 1 should have no prev link")

    big = get(client, f"{API_V1_PREFIX}/genes?per_page=50&page=1").json
    expect([r["gene_id"] for r in big["data"]][:25] == ids1,
           "ordering is not stable across different page sizes")


@case("documented example responses still match live responses", "contract")
def _examples_current(client):
    """A stale example is worse than none: a consumer will code against it.

    Regenerate with: python scripts/capture_api_examples.py
    """
    from capture_api_examples import capture_url

    problems = []
    for ep in ENDPOINTS:
        example = EXAMPLE_RESPONSES.get(ep["id"])
        if not example:
            continue
        r = get(client, capture_url(ep))
        if r.status != 200:
            problems.append(f"{ep['id']}: example URL now returns {r.status}")
            continue
        live = r.json

        missing_top = set(live) - set(example)
        if missing_top:
            problems.append(f"{ep['id']}: response gained top-level {sorted(missing_top)}")

        live_data, ex_data = live.get("data"), example.get("data")
        if isinstance(live_data, dict) and isinstance(ex_data, dict):
            missing = set(live_data) - set(ex_data)
            if missing:
                problems.append(f"{ep['id']}: data gained {sorted(missing)}")
        elif isinstance(live_data, list) and isinstance(ex_data, list):
            if live_data and ex_data:
                missing = set(live_data[0]) - set(ex_data[0])
                if missing:
                    problems.append(f"{ep['id']}: rows gained {sorted(missing)}")
            elif live_data and not ex_data:
                problems.append(f"{ep['id']}: example shows an empty list, live returns rows")

    expect(not problems,
           "example responses are stale — run scripts/capture_api_examples.py:\n    "
           + "\n    ".join(problems))


@case("the index is self-describing enough to use without other docs", "contract")
def _index_completeness(client):
    body = get(client, f"{API_V1_PREFIX}/").json
    for key in ("api_version", "base_url", "description", "conventions", "limits",
                "enums", "endpoints"):
        expect(key in body, f"index is missing {key!r}")
    documented = {e["id"] for e in ENDPOINTS if e["id"] != "index"}
    with_examples = {e["id"] for e in body["endpoints"] if "example_response" in e}
    expect(documented <= with_examples,
           f"endpoints listed without an example response: {sorted(documented - with_examples)}")
    for e in body["endpoints"]:
        expect(e["url_template"].startswith(body["base_url"]),
               f"{e['id']}: url_template {e['url_template']!r} is not under base_url")
        expect(e["example_url"].startswith("http"),
               f"{e['id']}: example_url {e['example_url']!r} is not absolute")
        expect(e["summary"], f"{e['id']}: no summary")
    compact = get(client, f"{API_V1_PREFIX}/?examples=false").json
    expect(all("example_response" not in e for e in compact["endpoints"]),
           "examples=false still returned example bodies")


DEFAULT_BODY_BUDGET = 16 * 1024
LINE_BUDGET = 5000


@case("default responses stay small and line-broken", "contract")
def _default_size(client):
    """Guards the regression that made /gene and /tf responses 80-230 KB: bulk
    lists shipped by default, serialised on a single line."""
    problems = []
    for ep in ENDPOINTS:
        if ep["id"] == "index":
            continue
        r = get(client, ep["example_url"])
        size = len(r.text.encode())
        longest = max(len(line) for line in r.text.split("\n"))
        if size > DEFAULT_BODY_BUDGET:
            problems.append(f"{ep['example_url']} is {size / 1024:.1f} KiB "
                            f"(budget {DEFAULT_BODY_BUDGET // 1024} KiB)")
        if longest > LINE_BUDGET:
            problems.append(f"{ep['example_url']} has a {longest}-character line")
    expect(not problems, "oversized default responses:\n    " + "\n    ".join(problems))

    # A table row is one line. Expanding every field onto its own line is what
    # inflated opt-in bulk responses to thousands of lines.
    r = get(client, f"{API_V1_PREFIX}/genes?per_page=50")
    rows = len(r.json["data"])
    lines = r.text.count("\n")
    expect(lines <= rows + 40,
           f"/genes?per_page=50 spans {lines} lines for {rows} rows — rows are being expanded")


@case("opt-in lists are counted when omitted and match when included", "contract")
def _counted_opt_ins(client):
    checks = [
        ("/gene/SOX17", "go_terms", "n_go_terms"),
        ("/gene/SOX17", "perturbation_effects", "n_perturbation_effects"),
        ("/gene/SOX17", "elements", "n_elements"),
        ("/module/DE-1.1", "genes", "n_genes"),
        ("/module/DE-1.1", "enrichment", "n_enrichment_terms"),
        ("/go-term/GO:0030183", "genes", "n_genes"),
    ]
    for path, field, count in checks:
        base = get(client, API_V1_PREFIX + path).json["data"]
        expect(field not in base, f"{path}: {field} should be opt-in")
        expect(isinstance(base.get(count), int), f"{path}: default response lacks {count}")
        full = get(client, f"{API_V1_PREFIX}{path}?include={field}").json["data"]
        expect(field in full, f"{path}?include={field} did not return {field}")
        truncated = full.get(f"{field}_truncated")
        limit = full.get(f"{field}_limit", full.get("genes_limit"))
        if truncated or (limit and base[count] > limit):
            expect(len(full[field]) <= base[count], f"{path}: {field} longer than {count}")
        else:
            expect(len(full[field]) == base[count],
                   f"{path}: {count}={base[count]} but include={field} returned "
                   f"{len(full[field])} rows")


@case("tf module regulation defaults to one collection and counts the rest", "contract")
def _regulation_levels(client):
    default = get(client, f"{API_V1_PREFIX}/tf/ARID1A").json["data"]
    expect(default["module_regulation_level"] == "hotspot_supermodule", "unexpected default level")
    expect(all(r["module_collection"] == "hotspot_supermodule"
               for r in default["module_regulation"]),
           "default level returned rows from other collections")
    counts = default["n_module_regulation"]
    expect(len(default["module_regulation"]) == counts["hotspot_supermodule"],
           "row count disagrees with n_module_regulation")
    everything = get(client, f"{API_V1_PREFIX}/tf/ARID1A?level=all").json["data"]
    expect(len(everything["module_regulation"]) == min(sum(counts.values()), 500),
           f"level=all returned {len(everything['module_regulation'])}, counts sum to "
           f"{sum(counts.values())}")


def _collect_urls(value):
    if isinstance(value, dict):
        for k, v in value.items():
            if k == "url_template":
                continue  # a template, not a link — /search needs q= filled in
            yield from _collect_urls(v)
    elif isinstance(value, list):
        for v in value:
            yield from _collect_urls(v)
    elif isinstance(value, str) and "/endoderm-perturbseq" in value:
        yield value


@case("every URL in a response is absolute and resolves", "contract")
def _absolute_links(client):
    """Chat assistants only fetch URLs they have already seen verbatim, so a
    relative or broken link is a dead end for them. Follow a sample to prove the
    API is navigable from its own responses."""
    from urllib.parse import urlsplit

    sources = ["/", "/gene/SOX17", "/tf/ARID1A", "/module/DE-1", "/genes?per_page=3",
               "/tfs?per_page=3", "/modules?per_page=3", "/search?q=SOX&limit=8",
               "/search?q=endoderm&type=go_term&limit=3",
               "/atac-peak/58938", "/atac-peak/58938/tfs?per_page=3", "/edges?per_page=3",
               "/link/tf-gene/FOXA2/SOX17", "/link/gene-module?module=DE-1&per_page=3",
               "/go-term/GO:0030183"]
    seen = set()
    for src in sources:
        body = get(client, API_V1_PREFIX + src).json
        for url in _collect_urls(body):
            expect(url.startswith(("http://", "https://")),
                   f"{src} contains a relative URL: {url!r}")
            if "{" not in url:
                seen.add(url)

    broken = []
    for url in sorted(seen):
        parts = urlsplit(url)
        path = parts.path + (f"?{parts.query}" if parts.query else "")
        r = client.request(path)
        if r.status != 200:
            broken.append(f"{r.status} {url}")
    expect(not broken, f"{len(broken)} of {len(seen)} links do not resolve:\n    "
           + "\n    ".join(broken[:15]))


@case("a forged Host header cannot redirect links elsewhere", "framework")
def _host_header(client):
    if not getattr(client, "app", None):
        return
    flask_client = client.app.test_client()
    r = flask_client.get(f"{API_V1_PREFIX}/gene/SOX17", base_url="https://evil.example")
    links = json.loads(r.get_data(as_text=True))["links"]
    expect(all(v.startswith("https://www.huangfulab.com/") for v in links.values() if v),
           f"forged Host leaked into links: {links}")
    r = flask_client.get(f"{API_V1_PREFIX}/gene/SOX17", base_url="https://huangfulab.com")
    links = json.loads(r.get_data(as_text=True))["links"]
    expect(links["self"].startswith("https://huangfulab.com/"),
           "a genuine host should be echoed so links match what the caller fetched")


@case("the JSON index describes every documented endpoint", "contract")
def _index(client):
    body = get(client, f"{API_V1_PREFIX}/").json
    expect(body["api_version"] == "v1", f"api_version was {body['api_version']!r}")
    expect(body["base_url"].startswith("http") and body["base_url"].endswith(API_V1_PREFIX),
           f"base_url should be an absolute URL, got {body['base_url']!r}")
    listed = {e["id"] for e in body["endpoints"]}
    expect(listed == {e["id"] for e in ENDPOINTS},
           f"index lists {listed}, catalog has {{e['id'] for e in ENDPOINTS}}")
    expect(body["limits"]["max_per_page"] == 500, "max_per_page not reported")


# ── structural cases (in-process only) ───────────────────────────────────────

@case("catalog covers exactly the registered routes", "framework")
def _coverage(client):
    if not getattr(client, "app", None):
        return
    registered = {
        r.endpoint for r in client.app.url_map.iter_rules()
        if r.endpoint.startswith("perturbseq_api_v1.")
    } - {"perturbseq_api_v1.unknown_endpoint"}
    documented = {e["endpoint"] for e in ENDPOINTS}
    missing = registered - documented
    stale = documented - registered
    expect(not missing, f"routes with no catalog entry: {sorted(missing)}")
    expect(not stale, f"catalog entries with no route: {sorted(stale)}")


@case("repeated requests do not leak file descriptors", "framework")
def _no_fd_leak(client):
    if not getattr(client, "app", None) or not os.path.isdir("/proc/self/fd"):
        return
    url = f"{API_V1_PREFIX}/gene/SOX17"
    get(client, url)
    before = len(os.listdir("/proc/self/fd"))
    for _ in range(200):
        get(client, url)
    delta = len(os.listdir("/proc/self/fd")) - before
    expect(delta < 5,
           f"file descriptors grew by {delta} over 200 requests — "
           f"api_v1_bp.teardown_request(close_db) is probably missing")


# Joining one of these without pinning the order risks a scan of 50M-185M rows.
_HUGE_TABLES = ("tf_gene_binding_scores", "atac_tf_overlaps", "tf_peaks", "gene_coexpression")
_UNPINNED_JOIN = re.compile(
    r"(?<!CROSS )\bJOIN\s+(" + "|".join(_HUGE_TABLES) + r")\b", re.IGNORECASE)


@case("joins against the 100M-row tables pin their order", "framework")
def _pinned_join_order(client):
    """A plan regression here is invisible in the response body, so it has to be
    checked in the source and then confirmed against the real planner.

    ANALYZE has never been run on this database, so SQLite plans these joins on
    default cardinality guesses. For tf_gene_binding_scores (184M rows) entering
    through the (dataset_id, gene_id) primary key takes ~0.008 s, while entering
    through idx_tf_bs_gene_id takes ~4.95 s for identical output.
    """
    problems = []
    for path in sorted((ROOT / "blueprints" / "api_v1").glob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
                continue
            if "plan-ok" in node.value:
                continue
            for table in _UNPINNED_JOIN.findall(node.value):
                problems.append(
                    f"{path.name}:{node.lineno} joins {table} without CROSS JOIN — "
                    f"the planner has no statistics and may pick a full scan")
    expect(not problems, "unpinned join order:\n    " + "\n    ".join(problems))

    if not getattr(client, "app", None):
        return
    import sqlite3

    from blueprints.perturbseq_bp import DB_PATH

    sql = (
        "SELECT bs.binding_score_A FROM tf_dataset_gene_table tdg "
        "CROSS JOIN tf_dataset_table td ON td.dataset_id = tdg.dataset_id "
        "CROSS JOIN tf_gene_binding_scores bs "
        "        ON bs.dataset_id = td.dataset_id AND bs.gene_id = ? "
        "WHERE tdg.gene_id = ? LIMIT 1"
    )
    db = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    try:
        plan = " | ".join(r[3] for r in db.execute("EXPLAIN QUERY PLAN " + sql, (1, 1)))
    finally:
        db.close()
    expect("idx_tf_bs_gene_id" not in plan,
           f"planner fell back to the gene-only index on a 184M-row table. Plan: {plan}")


_PAGED_SQL = re.compile(r"ORDER BY\s+(.+?)\s+LIMIT \? OFFSET \?", re.IGNORECASE | re.DOTALL)


@case("every paginated query orders by a unique column", "framework")
def _stable_paging(client):
    """Ordering by a non-unique column lets rows shift between pages, silently
    duplicating and dropping records. SQLite usually happens to be stable, so
    this cannot be detected from responses — it has to be checked in the source.
    """
    problems = []
    for path in sorted((ROOT / "blueprints" / "api_v1").glob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
                continue
            for order_by in _PAGED_SQL.findall(node.value):
                # A composite key can be unique without any single column being
                # so; /* stable-order */ marks one that was checked by hand.
                if "stable-order" in order_by:
                    continue
                last = order_by.split(",")[-1].strip().split()[0]
                if not last.endswith("_id"):
                    problems.append(
                        f"{path.name}:{node.lineno} pages on 'ORDER BY {order_by.strip()}' — "
                        f"the last term {last!r} is not unique. End on an _id column, or add "
                        f"/* stable-order */ if the composite really is total.")
    expect(not problems, "unstable pagination:\n    " + "\n    ".join(problems))


@case("an unmatched path under the prefix returns JSON, not HTML", "framework")
def _catch_all(client):
    u = f"{API_V1_PREFIX}/no/such/endpoint"
    body = assert_json_error(get(client, u), u, 404)
    expect(body["code"] == "unknown_endpoint", f"code was {body['code']!r}")


@case("the docs page carries no CORS header", "framework")
def _docs_not_cors(client):
    r = get(client, API_V1_PREFIX.rsplit("/", 1)[0])
    expect(r.headers.get("Access-Control-Allow-Origin") is None,
           "the HTML docs page should not advertise CORS")


_SQL_WORDS = re.compile(r"\b(SELECT|FROM|WHERE|JOIN|LIMIT|GROUP BY|ORDER BY)\b")
_WRITE_SQL = re.compile(r"\b(INSERT|UPDATE|DELETE|DROP|ALTER|ATTACH|CREATE|REPLACE)\s")
# Names permitted inside a SQL string. 'placeholders'/'ph' are the ','.join('?' * n)
# idiom; anything ALL_CAPS is a module-level SQL fragment constant, which by
# convention holds only literals. A lowercase local is rejected, because that is
# what a request-time value looks like.
_SAFE_INTERPOLATIONS = {"placeholders", "ph"}


def _is_safe_sql_name(name):
    if name is None:
        return False
    if name in _SAFE_INTERPOLATIONS:
        return True
    stripped = name.lstrip("_")
    return bool(stripped) and stripped.isupper()


def _root_name(node):
    while isinstance(node, (ast.Attribute, ast.Call, ast.Subscript)):
        node = node.func if isinstance(node, ast.Call) else node.value
    return node.id if isinstance(node, ast.Name) else None


@case("package source contains no write SQL, eval, jsonify or unsafe interpolation", "framework")
def _source_lint(client):
    problems = []
    for path in sorted((ROOT / "blueprints" / "api_v1").glob("*.py")):
        src = path.read_text()
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and \
                    node.func.id in ("eval", "exec"):
                problems.append(f"{path.name}:{node.lineno} calls {node.func.id}()")
            if isinstance(node, ast.Call) and _root_name(node.func) == "jsonify":
                problems.append(f"{path.name}:{node.lineno} calls jsonify() "
                                f"(NaN-unsafe; use json_response)")
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                if _WRITE_SQL.search(node.value) and _SQL_WORDS.search(node.value):
                    problems.append(f"{path.name}:{node.lineno} string contains write SQL")
            if isinstance(node, ast.JoinedStr):
                literal = "".join(v.value for v in node.values
                                  if isinstance(v, ast.Constant) and isinstance(v.value, str))
                if not _SQL_WORDS.search(literal):
                    continue
                if "sql-ok" in src.splitlines()[node.lineno - 1]:
                    continue
                for part in node.values:
                    if not isinstance(part, ast.FormattedValue):
                        continue
                    name = _root_name(part.value)
                    if not _is_safe_sql_name(name):
                        problems.append(
                            f"{path.name}:{node.lineno} f-string interpolates {name!r} "
                            f"into SQL")
            # Concatenation is the other way a value reaches SQL text, and an
            # f-string-only check would miss it entirely.
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
                operands = [node.left, node.right]
                literals = [o.value for o in operands
                            if isinstance(o, ast.Constant) and isinstance(o.value, str)]
                if not literals or not any(_SQL_WORDS.search(t) for t in literals):
                    continue
                if "sql-ok" in src.splitlines()[node.lineno - 1]:
                    continue
                for o in operands:
                    if isinstance(o, ast.Constant):
                        continue
                    name = _root_name(o) if not isinstance(o, ast.BinOp) else None
                    if isinstance(o, ast.BinOp):
                        continue
                    if not _is_safe_sql_name(name):
                        problems.append(
                            f"{path.name}:{node.lineno} concatenates {name!r} into SQL")
    expect(not problems, "source lint:\n    " + "\n    ".join(problems))


# ── runner ───────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-url", help="test a running server instead of an in-process app")
    ap.add_argument("--group", help="run only this group")
    ap.add_argument("--only", help="run only cases whose name contains this substring")
    ap.add_argument("--verbose", action="store_true", help="print the full assertion on failure")
    ap.add_argument("--list", action="store_true", help="list cases and exit")
    args = ap.parse_args()

    register_generated()

    selected = [
        (g, n, f) for g, n, f in CASES
        if (not args.group or g == args.group) and (not args.only or args.only in n)
    ]
    if args.list:
        for g, n, _ in selected:
            print(f"{g:12s} {n}")
        return 0
    if not selected:
        print("no cases matched")
        return 2

    if args.base_url:
        client = RemoteClient(args.base_url)
        target = args.base_url
    else:
        if not Path(DB_PATH).exists():
            print(f"DB not found at {DB_PATH}\nThis suite needs the real database.", file=sys.stderr)
            return 2
        client = LocalClient()
        target = "in-process"

    n_endpoints = len(ENDPOINTS)
    print(f"api/v1 test suite — {n_endpoints} endpoints, {len(selected)} cases ({target})\n")

    failures = []
    t0 = time.time()
    by_group = {}
    for group, name, fn in selected:
        by_group.setdefault(group, [])
        try:
            fn(client)
            by_group[group].append(".")
        except Exception as exc:
            by_group[group].append("F")
            failures.append((group, name, exc))

    for group, marks in by_group.items():
        ok = marks.count(".")
        bad = marks.count("F")
        tail = f"   {bad} FAILED" if bad else ""
        print(f"  {group:12s} {''.join(marks):<40s} {ok:3d} ok{tail}")

    if failures:
        print()
        for group, name, exc in failures:
            print(f"FAILED  {group} / {name}")
            text = str(exc) if args.verbose else str(exc).split("\n")[0][:300]
            print(f"  {text}")
            print()

    total = len(selected)
    print(f"\n{total} cases, {total - len(failures)} passed, {len(failures)} failed "
          f"in {time.time() - t0:.1f}s")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
