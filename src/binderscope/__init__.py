"""binderscope — structure-aware triage of protein binder designs.

A post-processing pipeline for binder design campaigns (BindCraft and
compatible outputs): restore an engineered target to wild type, score the
interface with PyRosetta, check where the binder actually sits relative to a
reference structure, rank designs on explicit weighted criteria, and render the
result as a self-contained dashboard.

Everything target-specific lives in a YAML config file.
"""

from .config import Config, load_config
from .mapping import SegmentMap
from .ranking import apply_filters, normalise, rank_designs

__version__ = "0.1.0"

__all__ = [
    "Config",
    "SegmentMap",
    "__version__",
    "apply_filters",
    "load_config",
    "normalise",
    "rank_designs",
]
