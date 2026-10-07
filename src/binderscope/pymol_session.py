"""Optional PyMOL session generation, one per design.

A metric table tells you which designs score well. A session you can open tells
you why, and whether the pose is believable. Sessions are written only when a
``pymol`` block is present in the config and the binary is found.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from .config import Config
from .mapping import SegmentMap

_SCRIPT = """\
load {complex_pdb}, {object_name}
load {reference_pdb}, reference
{align_command}
hide everything
show cartoon
bg_color white
color {binder_color}, {object_name} and chain {binder_chain}
color {target_color}, {object_name} and chain {target_chain}
{reference_colors}
dss
set seq_view, 1
save {session_path}
quit
"""


def find_pymol(config: Config) -> str | None:
    """Resolve the PyMOL executable from config or PATH."""
    configured = config.pymol.get("binary")
    if configured:
        path = Path(configured).expanduser()
        return str(path) if path.is_file() else None
    return shutil.which("pymol")


def write_session(
    complex_pdb: str | Path,
    session_path: str | Path,
    config: Config,
    segment_map: SegmentMap | None = None,
    pymol_binary: str | None = None,
    timeout: int = 180,
) -> bool:
    """Render one session. Returns False when PyMOL is unavailable or failed."""
    binary = pymol_binary or find_pymol(config)
    if not binary or not config.reference:
        return False

    complex_pdb = Path(complex_pdb)
    session_path = Path(session_path)
    session_path.parent.mkdir(parents=True, exist_ok=True)
    object_name = complex_pdb.stem.replace("-", "_").replace(".", "_")

    align_command = ""
    if segment_map and segment_map.segments:
        anchor = segment_map.segments[0]
        reference_end = anchor.reference_start + anchor.length - 1
        align_command = (
            f"align {object_name} and chain {config.target_chain} and "
            f"resi {anchor.input_start}-{anchor.input_end}, "
            f"reference and chain {config.reference.target_chain} and "
            f"resi {anchor.reference_start}-{reference_end}"
        )

    palette = config.pymol.get("colors", {})
    reference_colors = "\n".join(
        f"color {palette.get(target.label, 'grey70')}, "
        f"reference and chain {target.chain}"
        for target in config.reference.clash_targets
    )

    script = _SCRIPT.format(
        complex_pdb=complex_pdb,
        object_name=object_name,
        reference_pdb=config.reference.structure,
        align_command=align_command,
        binder_chain=config.binder_chain,
        target_chain=config.target_chain,
        binder_color=palette.get("binder", "cyan"),
        target_color=palette.get("target", "wheat"),
        reference_colors=reference_colors,
        session_path=session_path,
    )

    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as fh:
        fh.write(script)
        script_path = fh.name
    try:
        result = subprocess.run(
            [binary, "-cq", script_path],
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        return result.returncode == 0 and session_path.is_file()
    except (subprocess.TimeoutExpired, OSError):
        return False
    finally:
        Path(script_path).unlink(missing_ok=True)
