"""Public read-only JSON API for the Endoderm Perturb-Seq Browser.

Mounted at /endoderm-perturbseq/api/v1, with its human documentation at
/endoderm-perturbseq/api. Deliberately separate from the internal /api/* routes
in perturbseq_bp.py, which are page plumbing and free to change.
"""
from .core import api_docs_bp, api_v1_bp, is_api_v1_path, routing_error_response
from . import (  # noqa: F401  registers the views
    routes_objects,
    routes_links,
    routes_collections,
    meta,
)

__all__ = ["api_v1_bp", "api_docs_bp", "is_api_v1_path", "routing_error_response"]
