"""The API's own metadata: a JSON index, the HTML docs page, and the catch-all
that turns an unmatched path into a JSON 404 instead of an HTML one.

api_docs_bp does not touch the database. If that ever changes it needs its own
teardown_request, or it will leak connections.
"""
from flask import redirect, render_template

from .catalog import ENDPOINTS, ENUMS
from .examples import EXAMPLE_RESPONSES
from .core import (
    API_DOCS_PATH,
    API_V1_PREFIX,
    ApiError,
    MAX_RESPONSE_BYTES,
    PER_PAGE_DEFAULT,
    PER_PAGE_MAX,
    QUERY_DEADLINE_COLLECTION,
    QUERY_DEADLINE_OBJECT,
    _LAST_COMMIT_TS,
    _fmt_time_ago,
    absolute,
    api_docs_bp,
    api_v1_bp,
    arg_bool,
    json_response,
    public_origin,
    reject_unknown_args,
)

_LIMITS = {
    "default_per_page": PER_PAGE_DEFAULT,
    "max_per_page": PER_PAGE_MAX,
    "max_response_bytes": MAX_RESPONSE_BYTES,
    "query_timeout_seconds": {
        "objects_and_links": QUERY_DEADLINE_OBJECT,
        "collections_and_search": QUERY_DEADLINE_COLLECTION,
    },
}

_CONVENTIONS = {
    "entity_response":
        "A single object returns {\"data\": {...}, \"links\": {...}} and, when there is "
        "something to report, \"meta\". 'links' sits beside 'data', never inside it, so "
        "'data' flattens cleanly into a table or data frame.",
    "collection_response":
        "A listing returns {\"data\": [...], \"page\", \"per_page\", \"total\", \"pages\", "
        "\"total_is_exact\", \"links\"}. Follow links.next until it is null. A 'total' of "
        "null means the count was too expensive to compute, never that there are no rows.",
    "errors":
        "Always JSON: {\"error\": \"...\", \"status\": 404, \"code\": \"gene_not_found\", "
        "\"path\": \"...\"}. Branch on 'code', which is stable; 'error' is prose for humans "
        "and may be reworded.",
    "pagination":
        "1-based 'page' with 'per_page'. A per_page above the maximum is rejected with 400 "
        "rather than silently reduced, so a short page always means you reached the end.",
    "missing_and_infinite_values":
        "Absent measurements are null. Where a Fisher test separated perfectly the odds "
        "ratio is null with \"odds_ratio_infinite\": true beside it, because bare Infinity "
        "is not valid JSON and a sentinel number would be mistaken for a measurement.",
    "truncation":
        "Bounded lists carry a companion \"<field>_truncated\" boolean, so a capped list is "
        "always distinguishable from a genuinely short one.",
    "name_resolution":
        "Genes accept a symbol, an Ensembl gene ID or a synonym; when the input was not the "
        "canonical symbol, meta.resolved_from and meta.match_type say what matched. Gene "
        "clusters accept GC1, TC-1, gene_cluster_1 or cluster_1 and always answer with GC1.",
    "urls":
        "Every URL in a response is absolute and can be fetched as-is. Chat assistants "
        "that only fetch URLs already seen in the conversation can therefore follow links, "
        "row 'url' fields and pagination from any response they have read.",
    "optional_blocks":
        "Entity responses are kept small. Long lists are opt-in via include= and are "
        "represented in the default response by a count named n_<field>.",
    "example_responses":
        "Captured from the live API and trimmed: arrays are cut to a single element and "
        "long strings to ~140 characters. Field names and types are verbatim. Examples "
        "are captured with every include= block enabled, so they show the fullest shape.",
}

_GROUP_ORDER = ("objects", "links", "collections", "meta")
_GROUP_TITLES = {
    "objects": "Objects",
    "links": "Links",
    "collections": "Collections & search",
    "meta": "Meta",
}


def _index_entry(endpoint, examples):
    path = API_V1_PREFIX + (endpoint["path"] if endpoint["path"] != "/" else "/")
    entry = {
        "id": endpoint["id"],
        "group": endpoint["group"],
        "method": endpoint["method"],
        "url_template": absolute(path),
        "summary": endpoint["summary"],
        "path_params": endpoint["path_params"],
        "query_params": endpoint["query_params"],
        "returns": endpoint["returns"],
        "example_url": absolute(endpoint["example_url"]),
        "errors": endpoint["errors"],
        "notes": endpoint["notes"],
    }
    if examples and endpoint["id"] in EXAMPLE_RESPONSES:
        entry["example_response"] = EXAMPLE_RESPONSES[endpoint["id"]]
    return entry


@api_v1_bp.route("")
@api_v1_bp.route("/")
def api_index():
    """The whole API described in one document.

    Example responses are included by default: this URL is meant to be handed to
    a person or a program that knows nothing about the API, and a spec that omits
    response shapes invites a consumer to invent them.
    """
    reject_unknown_args("examples")
    examples = arg_bool("examples", True)
    return json_response({
        "api_version": "v1",
        "base_url": absolute(API_V1_PREFIX),
        "docs": absolute(API_DOCS_PATH),
        "data_updated": _fmt_time_ago(_LAST_COMMIT_TS),
        "description":
            "Read-only JSON API for the Endoderm Perturb-Seq Browser, a CRISPR "
            "transcription-factor perturbation screen in human definitive endoderm "
            "differentiation. No authentication and no rate limit. Every endpoint is GET. "
            "This document lists every endpoint, its parameters, and an example response.",
        "conventions": _CONVENTIONS,
        "limits": _LIMITS,
        "enums": {k: list(v) for k, v in ENUMS.items()},
        "endpoints": [_index_entry(e, examples) for e in ENDPOINTS],
    })


@api_v1_bp.route("/<path:unknown>")
def unknown_endpoint(unknown):
    """Blueprint error handlers never fire for unmatched URLs — on a routing
    failure request.blueprints is empty — so an explicit catch-all is the only
    way to return JSON here. Werkzeug ranks static rules above <path:>, so this
    shadows nothing."""
    if unknown.endswith("/"):
        return redirect(f"{API_V1_PREFIX}/{unknown.rstrip('/')}", code=308)
    raise ApiError(
        404,
        f"Unknown endpoint '/{unknown}'. See {absolute(API_V1_PREFIX)}/ for the endpoint index.",
        "unknown_endpoint",
    )


@api_docs_bp.route("")
@api_docs_bp.route("/")
def api_docs():
    groups = [
        (g, _GROUP_TITLES[g], [e for e in ENDPOINTS if e["group"] == g])
        for g in _GROUP_ORDER
    ]
    # Absolute for the same reason as the JSON: an assistant reading this page
    # can only follow links that appear on it verbatim.
    return render_template(
        "perturbseq/api.html",
        groups=[grp for grp in groups if grp[2]],
        enums=ENUMS,
        origin=public_origin(),
        base=absolute(API_V1_PREFIX),
        limits=_LIMITS,
    )
