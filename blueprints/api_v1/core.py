"""Scaffolding for the public API: blueprints, response envelope, parameter
validation, error handling, CORS, and the read-only DB guard.

Nothing in this package calls jsonify — json_response() is the only place a body
is produced, and it uses allow_nan=False so a NaN or Infinity that escaped
rows_to_dicts raises here instead of reaching a consumer as invalid JSON.
"""
import json
import math
import re
import sqlite3
import time
from datetime import date, datetime
from decimal import Decimal
from urllib.parse import quote, urlencode

from flask import Blueprint, Response, current_app, g, request
from werkzeug.exceptions import HTTPException

from blueprints.perturbseq_bp import (
    PERTURBSEQ_PREFIX,
    MODULE_COLORS,
    _LAST_COMMIT_TS,
    _db_maintenance_active,
    _fmt_time_ago,
    close_db,
    get_db,
    rows_to_dicts,
)
from .catalog import ENUMS, P_PER_PAGE

API_V1_PREFIX = PERTURBSEQ_PREFIX + "/api/v1"
API_DOCS_PATH = PERTURBSEQ_PREFIX + "/api"

PER_PAGE_MAX = 500
# Deliberately small: agent fetch tools refuse or truncate large bodies, and a
# default page is what an unfamiliar caller sees first. Scripts that want bulk
# data raise per_page themselves.
PER_PAGE_DEFAULT = 25
PAGE_MAX = 100_000
MAX_RESPONSE_BYTES = 5 * 1024 * 1024
QUERY_DEADLINE_OBJECT = 5.0
QUERY_DEADLINE_COLLECTION = 2.0

assert P_PER_PAGE["max"] == PER_PAGE_MAX, "catalog and core disagree on the page-size cap"
assert P_PER_PAGE["default"] == PER_PAGE_DEFAULT, "catalog and core disagree on the page size"

# api_v1_bp carries the JSON contract: CORS, JSON errors, the read-only guard.
# api_docs_bp serves one HTML page and deliberately shares none of that.
api_v1_bp = Blueprint("perturbseq_api_v1", __name__, url_prefix=API_V1_PREFIX)
api_docs_bp = Blueprint("perturbseq_api_docs", __name__, url_prefix=API_DOCS_PATH)

# close_db is registered on perturbseq_bp, whose teardown does not fire for this
# blueprint. Without this line every API request leaks a handle to a 33 GB database.
api_v1_bp.teardown_request(close_db)


class ApiError(Exception):
    def __init__(self, status, message, code="invalid_param"):
        super().__init__(message)
        self.status = status
        self.message = message
        self.code = code


# ── serialization ────────────────────────────────────────────────────────────

def _json_default(o):
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    if isinstance(o, (date, datetime)):
        return o.isoformat()
    if isinstance(o, Decimal):
        return float(o)
    if hasattr(o, "item"):
        return o.item()
    raise TypeError(f"unserialisable {type(o).__name__}")


def scrub(value):
    """NaN/Inf sanitisation for scalars computed in Python.

    Rows out of SQLite go through rows_to_dicts, which applies the same rules.
    """
    if isinstance(value, float):
        if math.isnan(value):
            return None
        if math.isinf(value):
            return 1e308 if value > 0 else -1e308
    return value


_LINE_WIDTH = 100
_FLAT_LIST_WIDTH = 400


def _format(value, level):
    """Pretty JSON that keeps records on one line.

    Plain indent=2 inflates payloads by ~50% (every scalar in every row gets its
    own line), which pushed responses past what agent fetch tools will read.
    Fully minified output has the opposite problem: a whole response on a single
    80,000-character line. Here an object of scalars is a record and always
    renders as one line — a table row per line — and anything else expands only
    when it does not fit within _LINE_WIDTH.
    """
    flat = json.dumps(value, ensure_ascii=False, separators=(", ", ": "))
    pad = "  " * level
    if not isinstance(value, (dict, list)) or not value or len(pad) + len(flat) <= _LINE_WIDTH:
        return flat
    children = value.values() if isinstance(value, dict) else value
    if not any(isinstance(c, (dict, list)) and c for c in children):
        if isinstance(value, dict) or len(flat) <= _FLAT_LIST_WIDTH:
            return flat
    inner = "  " * (level + 1)
    if isinstance(value, dict):
        items = [f"{inner}{json.dumps(k, ensure_ascii=False)}: {_format(v, level + 1)}"
                 for k, v in value.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    items = [inner + _format(v, level + 1) for v in value]
    return "[\n" + ",\n".join(items) + "\n" + pad + "]"


def json_response(payload, status=200):
    # Round-trip first: this applies _json_default and allow_nan=False once, so a
    # NaN or Infinity that escaped rows_to_dicts raises here rather than being
    # written into the formatted body.
    plain = json.loads(json.dumps(payload, allow_nan=False, default=_json_default))
    body = (_format(plain, 0) + "\n").encode()
    if status == 200 and len(body) > MAX_RESPONSE_BYTES:
        return error_response(
            413,
            f"Response would be {len(body)} bytes, over the {MAX_RESPONSE_BYTES} byte limit. "
            f"Lower per_page, narrow your filters, or drop an include= block.",
            "payload_too_large",
        )
    resp = Response(body, status=status, mimetype="application/json")
    resp.headers["Content-Length"] = str(len(body))
    return resp


def error_response(status, message, code):
    return json_response(
        {"error": message, "status": status, "code": code, "path": request.path},
        status=status,
    )


def _base_meta(extra=None):
    meta = dict(extra or {})
    if _db_maintenance_active():
        meta["db_maintenance"] = True
    return meta


def entity(data, links=None, meta=None):
    payload = {"data": data, "links": links or {}}
    meta = _base_meta(meta)
    if meta:
        payload["meta"] = meta
    return json_response(payload)


def collection(rows, *, page, per_page, total, total_is_exact=True, meta=None, paginated=True):
    if paginated:
        links = _page_links(page, per_page, payload_total=total, n_rows=len(rows))
    else:
        # An endpoint that does not accept page/per_page must not advertise
        # links carrying them: following one would return a 400.
        query_string = request.query_string.decode()
        self_url = absolute(request.path + (f"?{query_string}" if query_string else ""))
        links = {"self": self_url, "first": None, "prev": None, "next": None, "last": None}
    payload = {
        "data": rows,
        "page": page,
        "per_page": per_page,
        "total": total,
        "total_is_exact": total_is_exact,
        "pages": (math.ceil(total / per_page) if total else 0) if total is not None else None,
        "links": links,
    }
    meta = _base_meta(meta)
    if total is not None and not total_is_exact:
        meta["count_note"] = f"Exact count exceeds {total}; 'total' is a lower bound."
    if meta:
        payload["meta"] = meta
    return json_response(payload)


def _page_url(page, per_page):
    args = {k: v for k, v in request.args.items() if k not in ("page", "per_page")}
    args["page"] = page
    args["per_page"] = per_page
    return f"{public_origin()}{request.path}?{urlencode(sorted(args.items()))}"


def _page_links(page, per_page, *, payload_total, n_rows):
    last = math.ceil(payload_total / per_page) if payload_total else None
    has_next = (page < last) if last is not None else (n_rows == per_page)
    return {
        "self": _page_url(page, per_page),
        "first": _page_url(1, per_page),
        "prev": _page_url(page - 1, per_page) if page > 1 else None,
        "next": _page_url(page + 1, per_page) if has_next else None,
        "last": _page_url(last, per_page) if last else None,
    }


CANONICAL_ORIGIN = "https://www.huangfulab.com"
_PUBLIC_HOSTS = frozenset({"huangfulab.com", "www.huangfulab.com"})
_LOCAL_HOSTS = frozenset({"localhost", "127.0.0.1"})


def public_origin():
    """Scheme and host for the absolute URLs in every response.

    URLs are absolute because chat assistants such as claude.ai will only fetch a
    URL that has already appeared verbatim in the conversation or an earlier
    result. A relative path never qualifies, so an agent that has read one
    response could not follow any of its links.

    The request's own host is echoed only when it is one of ours, so a link keeps
    the exact form (with or without www) that the caller fetched. Any other Host
    header falls back to the canonical origin: echoing it would let a forged
    header plant links to another site in a publicly cacheable response.
    """
    host = (request.host or "").lower()
    name = host.split(":", 1)[0]
    if name in _PUBLIC_HOSTS:
        return f"https://{name}"
    if name in _LOCAL_HOSTS:
        return f"{request.scheme}://{host}"
    return CANONICAL_ORIGIN


def absolute(path):
    return public_origin() + path


def path_link(*segments):
    return absolute(API_V1_PREFIX + "".join("/" + quote(str(s), safe="") for s in segments))


def html_link(*segments):
    return absolute(PERTURBSEQ_PREFIX + "".join("/" + quote(str(s), safe="") for s in segments))


# ── parameter validation ─────────────────────────────────────────────────────

def arg_int(name, default=None, *, minimum=None, maximum=None):
    if name not in request.args:
        return default
    # type=int returns None on garbage rather than raising; the membership test
    # above is what separates "absent" from "malformed".
    val = request.args.get(name, type=int)
    if val is None:
        raise ApiError(400, f"Parameter '{name}' must be an integer, got {request.args[name]!r}.")
    if minimum is not None and val < minimum:
        raise ApiError(400, f"Parameter '{name}' must be >= {minimum}, got {val}.")
    if maximum is not None and val > maximum:
        raise ApiError(400, f"Parameter '{name}' must be <= {maximum}, got {val}.")
    return val


def arg_float(name, default=None, *, minimum=None, maximum=None):
    if name not in request.args:
        return default
    val = request.args.get(name, type=float)
    if val is None or val != val or math.isinf(val):
        raise ApiError(400, f"Parameter '{name}' must be a finite number, "
                            f"got {request.args[name]!r}.")
    if minimum is not None and val < minimum:
        raise ApiError(400, f"Parameter '{name}' must be >= {minimum}, got {val}.")
    if maximum is not None and val > maximum:
        raise ApiError(400, f"Parameter '{name}' must be <= {maximum}, got {val}.")
    return val


def arg_str(name, default="", *, max_len=120):
    val = request.args.get(name, default)
    if val is None:
        return default
    val = val.strip()
    if len(val) > max_len:
        raise ApiError(400, f"Parameter '{name}' must be at most {max_len} characters.")
    return val


_TRUE = {"true", "1", "yes", "on"}
_FALSE = {"false", "0", "no", "off"}


def arg_bool(name, default=None):
    if name not in request.args:
        return default
    raw = request.args[name].strip().lower()
    if raw in _TRUE:
        return True
    if raw in _FALSE:
        return False
    raise ApiError(400, f"Parameter '{name}' must be true or false, got {request.args[name]!r}.")


def arg_enum(name, allowed, default=None):
    if name not in request.args:
        return default
    val = request.args[name].strip()
    if val not in allowed:
        raise ApiError(400, f"Parameter '{name}' must be one of {', '.join(allowed)}; got {val!r}.")
    return val


def arg_enum_list(name, allowed):
    raw = request.args.get(name, "").strip()
    if not raw:
        return set()
    values = {v.strip() for v in raw.split(",") if v.strip()}
    unknown = values - set(allowed)
    if unknown:
        raise ApiError(
            400,
            f"Parameter '{name}' has unknown value(s) {', '.join(sorted(unknown))}; "
            f"valid values are {', '.join(allowed)}.",
        )
    return values


_NAME_OK = re.compile(r"^[A-Za-z0-9:._\-+ ()/]{1,120}$")


def clean_name(value, what):
    value = (value or "").strip()
    if not value:
        raise ApiError(400, f"A {what} name is required.")
    if not _NAME_OK.match(value):
        raise ApiError(400, f"{what.capitalize()} name {value!r} contains unsupported characters.")
    return value


def page_params():
    page = arg_int("page", 1, minimum=1, maximum=PAGE_MAX)
    per_page = arg_int("per_page", PER_PAGE_DEFAULT, minimum=1, maximum=PER_PAGE_MAX)
    return page, per_page, (page - 1) * per_page


def reject_unknown_args(*known):
    unknown = set(request.args) - set(known)
    if unknown:
        raise ApiError(
            400,
            f"Unknown query parameter(s): {', '.join(sorted(unknown))}. "
            f"This endpoint accepts: {', '.join(sorted(known)) or 'none'}.",
            "unknown_param",
        )


# ── database access ──────────────────────────────────────────────────────────

_ALLOWED_ACTIONS = frozenset({
    sqlite3.SQLITE_SELECT,
    sqlite3.SQLITE_READ,
    sqlite3.SQLITE_FUNCTION,
})


def _readonly_authorizer(action, arg1, arg2, db_name, trigger):
    return sqlite3.SQLITE_OK if action in _ALLOWED_ACTIONS else sqlite3.SQLITE_DENY


def api_db(deadline=QUERY_DEADLINE_OBJECT):
    """Per-request connection, structurally restricted to reads and time-boxed.

    The authorizer is stronger than a mode=ro connection string — it also blocks
    ATTACH and PRAGMA — and, unlike mode=ro, it cannot fail on a WAL database
    whose -shm file the process may not be able to write.
    """
    db = get_db()
    if not g.get("_api_v1_guarded"):
        db.set_authorizer(_readonly_authorizer)
        g._api_v1_guarded = True
    end = time.monotonic() + deadline
    db.set_progress_handler(lambda: 1 if time.monotonic() > end else 0, 10_000)
    return db


def query(db, sql, params=()):
    return rows_to_dicts(db.execute(sql, params).fetchall())


def query_one(db, sql, params=()):
    rows = query(db, sql, params)
    return rows[0] if rows else None


def scalar(db, sql, params=(), default=None):
    row = db.execute(sql, params).fetchone()
    return row[0] if row is not None and row[0] is not None else default


def _is_interrupt(exc):
    return isinstance(exc, sqlite3.OperationalError) and "interrupt" in str(exc).lower()


# ── request lifecycle ────────────────────────────────────────────────────────

@api_v1_bp.before_request
def _maintenance_gate():
    if _db_maintenance_active():
        resp = error_response(
            503,
            "The Perturb-Seq database is being updated and the API is temporarily offline. "
            "Try again shortly.",
            "db_maintenance",
        )
        resp.headers["Retry-After"] = "900"
        return resp
    return None


@api_v1_bp.after_request
def _api_headers(resp):
    resp.headers["Access-Control-Allow-Origin"] = "*"
    resp.headers["Access-Control-Allow-Methods"] = "GET, HEAD, OPTIONS"
    resp.headers["Access-Control-Allow-Headers"] = "Content-Type, Accept"
    resp.headers["Access-Control-Expose-Headers"] = "Content-Length, Retry-After"
    resp.headers["Access-Control-Max-Age"] = "86400"
    resp.headers["X-Content-Type-Options"] = "nosniff"
    resp.headers.setdefault("Cache-Control", "public, max-age=300")
    return resp


_HTTP_MESSAGES = {
    400: ("Bad request.", "bad_request"),
    404: ("Not found.", "not_found"),
    405: ("Method not allowed. All API endpoints are GET.", "method_not_allowed"),
    413: ("Response too large.", "payload_too_large"),
    503: ("The Perturb-Seq database is temporarily unavailable. Try again shortly.",
          "db_unavailable"),
}


@api_v1_bp.errorhandler(ApiError)
def _handle_api_error(e):
    return error_response(e.status, e.message, e.code)


@api_v1_bp.errorhandler(HTTPException)
def _handle_http_error(e):
    message, code = _HTTP_MESSAGES.get(e.code, ("Request could not be completed.", "http_error"))
    return error_response(e.code or 500, message, code)


@api_v1_bp.errorhandler(Exception)
def _handle_unexpected(e):
    if _is_interrupt(e):
        return error_response(
            503,
            "Query exceeded its time limit. Narrow your filters, lower per_page, "
            "or drop an include= block.",
            "query_timeout",
        )
    current_app.logger.exception("api/v1 failed: %s", request.path)
    # Deliberately a constant: str(e) here would leak SQL and paths.
    return error_response(500, "Internal server error.", "internal_error")


def is_api_v1_path():
    return request.path == API_V1_PREFIX or request.path.startswith(API_V1_PREFIX + "/")


def routing_error_response(status, message, code):
    """Used by app.py for 405, which Werkzeug raises before any blueprint is
    resolved, so no blueprint-scoped handler can see it."""
    resp = error_response(status, message, code)
    resp.headers["Access-Control-Allow-Origin"] = "*"
    return resp


@api_docs_bp.context_processor
def _docs_context():
    """base.html reads these from perturbseq_bp's context processor, which is
    blueprint-scoped. Without this mirror the maintenance banner silently
    vanishes on the docs page — Jinja renders the Undefined as falsy."""
    return {
        "module_colors": MODULE_COLORS,
        "last_updated": _fmt_time_ago(_LAST_COMMIT_TS),
        "db_maintenance": _db_maintenance_active(),
    }
