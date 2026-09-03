#!/usr/bin/env python3
"""Apply/inspect/restore the CRD_DAP cached-whitener patch for pPXF 9.4.8.

This utility deliberately edits only the *installed* ``ppxf/ppxf.py`` file in
whatever Python/conda environment is executing this script.  It does not ship a
copy of pPXF source code.

The patch adds one optional pPXF keyword::

    noise_inv_cholesky=None

When provided, it must be the precomputed inverse lower-Cholesky factor W=L^-1
of a covariance matrix C=L L^T.  pPXF is still passed a valid ordinary 1-D
``noise`` vector for its normal API validation, but after that validation the
patch installs W directly as pPXF's internal whitening operator.  This skips
repeating the covariance Cholesky factorization for every fit that reuses the
same covariance matrix.

The patch is intentionally version-locked to pPXF 9.4.8 and refuses to modify
an unexpected source layout.
"""
from __future__ import annotations

import argparse
import ast
import importlib
import inspect
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

PATCH_VERSION = "CRD_DAP_PRECOMPUTED_WHITENER_PATCH_V1"
SUPPORTED_PPXF_VERSION = "9.4.8"
BACKUP_SUFFIX = ".crd_dap_pre_cached_whitener_9_4_8.bak"


def _import_ppxf():
    try:
        import ppxf
        import ppxf.ppxf as ppxf_module
    except Exception as exc:
        raise RuntimeError(
            "Could not import pPXF from the active Python environment. "
            "Activate the CRD_DAP conda environment and install ppxf==9.4.8 first."
        ) from exc
    return ppxf, ppxf_module


def _ppxf_version(ppxf_pkg) -> str:
    value = getattr(ppxf_pkg, "__version__", None)
    if value is None:
        raise RuntimeError("Installed pPXF does not expose __version__; refusing to patch.")
    return str(value)


def _source_path(ppxf_module) -> Path:
    path = Path(inspect.getfile(ppxf_module)).expanduser().resolve()
    if path.name != "ppxf.py":
        raise RuntimeError(f"Expected installed module to be ppxf.py, got: {path}")
    return path


def _is_inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def _find_ppxf_init(tree: ast.AST) -> ast.FunctionDef:
    for node in getattr(tree, "body", []):
        if isinstance(node, ast.ClassDef) and node.name == "ppxf":
            for child in node.body:
                if isinstance(child, ast.FunctionDef) and child.name == "__init__":
                    return child
    raise RuntimeError("Could not find class ppxf.__init__ in installed ppxf.py.")


def _call_attr_names(node: ast.AST) -> set[str]:
    names: set[str] = set()
    for sub in ast.walk(node):
        if isinstance(sub, ast.Call):
            func = sub.func
            if isinstance(func, ast.Attribute):
                names.add(func.attr)
            elif isinstance(func, ast.Name):
                names.add(func.id)
    return names


def _find_covariance_preprocessing_if(init_node: ast.FunctionDef) -> ast.If:
    """Locate the top-level __init__ If that performs covariance factorization."""
    matches: list[ast.If] = []
    for stmt in init_node.body:
        if isinstance(stmt, ast.If):
            names = _call_attr_names(stmt)
            if "cholesky" in names and "solve_triangular" in names:
                matches.append(stmt)
    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one top-level covariance preprocessing block containing "
            f"cholesky + solve_triangular; found {len(matches)}. Source layout is not the "
            "validated pPXF-9.4.8 structure, so no modification was made."
        )
    return matches[0]


def _signature_has_keyword(init_node: ast.FunctionDef, name: str) -> bool:
    return name in [arg.arg for arg in init_node.args.args + init_node.args.kwonlyargs]


def _patched_source(source: str) -> str:
    if PATCH_VERSION in source:
        return source

    tree = ast.parse(source)
    init_node = _find_ppxf_init(tree)
    if _signature_has_keyword(init_node, "noise_inv_cholesky"):
        raise RuntimeError(
            "pPXF already has a noise_inv_cholesky argument but does not contain the CRD_DAP "
            "patch marker. Refusing to guess how this installation was modified."
        )

    cov_if = _find_covariance_preprocessing_if(init_node)
    if cov_if.end_lineno is None:
        raise RuntimeError("Python AST did not report end_lineno for covariance block.")

    lines = source.splitlines(keepends=True)

    # The function header spans from def-line through the line before the first
    # body statement. Add the keyword just before the final '):'. This avoids
    # depending on pPXF's line wrapping while keeping the rest of the file intact.
    first_body_line = min(stmt.lineno for stmt in init_node.body)
    header_start = init_node.lineno - 1
    header_end = first_body_line - 1
    header = "".join(lines[header_start:header_end])
    close = header.rfind("):")
    if close < 0:
        raise RuntimeError("Could not identify the end of ppxf.__init__ signature.")
    header = header[:close] + ", noise_inv_cholesky=None" + header[close:]
    lines[header_start:header_end] = [header]

    # Header replacement can change physical line count only if the source had
    # unusual newline handling. Re-parse to locate the covariance block again.
    interim = "".join(lines)
    tree2 = ast.parse(interim)
    init2 = _find_ppxf_init(tree2)
    cov_if2 = _find_covariance_preprocessing_if(init2)
    if cov_if2.end_lineno is None:
        raise RuntimeError("Could not re-locate covariance block after signature edit.")

    insert_at = cov_if2.end_lineno
    indent = " " * 8
    block = (
        f"\n{indent}# {PATCH_VERSION}\n"
        f"{indent}# Optional exact reuse of pPXF's already-computed whitening operator.\n"
        f"{indent}if noise_inv_cholesky is not None:\n"
        f"{indent}    noise_inv_cholesky = np.asarray(noise_inv_cholesky, dtype=float)\n"
        f"{indent}    if galaxy.ndim != 1:\n"
        f"{indent}        raise ValueError('noise_inv_cholesky is supported only for one-dimensional galaxy spectra')\n"
        f"{indent}    expected_shape = (galaxy.shape[0], galaxy.shape[0])\n"
        f"{indent}    if noise_inv_cholesky.shape != expected_shape:\n"
        f"{indent}        raise ValueError(f'noise_inv_cholesky must have shape {{expected_shape}}')\n"
        f"{indent}    if not np.all(np.isfinite(noise_inv_cholesky)):\n"
        f"{indent}        raise ValueError('noise_inv_cholesky must contain only finite values')\n"
        f"{indent}    if np.any(np.diag(noise_inv_cholesky) <= 0):\n"
        f"{indent}        raise ValueError('noise_inv_cholesky must have a strictly positive diagonal')\n"
        f"{indent}    if not np.allclose(noise_inv_cholesky, np.tril(noise_inv_cholesky), rtol=0.0, atol=1e-12):\n"
        f"{indent}        raise ValueError('noise_inv_cholesky must be the inverse of a lower-triangular Cholesky factor')\n"
        f"{indent}    self.noise = noise_inv_cholesky\n"
        f"{indent}    self.crd_dap_noise_inv_cholesky = True\n"
        f"{indent}else:\n"
        f"{indent}    self.crd_dap_noise_inv_cholesky = False\n"
    )

    lines2 = interim.splitlines(keepends=True)
    lines2.insert(insert_at, block)
    patched = "".join(lines2)

    compile(patched, "ppxf.py", "exec")
    tree3 = ast.parse(patched)
    init3 = _find_ppxf_init(tree3)
    if not _signature_has_keyword(init3, "noise_inv_cholesky"):
        raise RuntimeError("Internal patch verification failed: keyword missing from signature.")
    if PATCH_VERSION not in patched:
        raise RuntimeError("Internal patch verification failed: patch marker missing.")
    return patched


def _status(path: Path, source: str, version: str) -> None:
    print(f"Python executable : {sys.executable}")
    print(f"sys.prefix        : {Path(sys.prefix).resolve()}")
    print(f"pPXF version      : {version}")
    print(f"ppxf.py           : {path}")
    print(f"inside sys.prefix : {_is_inside(path, Path(sys.prefix))}")
    print(f"patch installed   : {PATCH_VERSION in source}")
    print(f"backup            : {Path(str(path) + BACKUP_SUFFIX)}")


def _subprocess_smoke() -> None:
    code = (
        "import inspect, ppxf; "
        "from ppxf.ppxf import ppxf as cls; "
        f"assert str(ppxf.__version__) == '{SUPPORTED_PPXF_VERSION}'; "
        "assert 'noise_inv_cholesky' in inspect.signature(cls).parameters"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Patch only the active environment's ppxf/ppxf.py to accept a cached inverse-Cholesky whitener."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="Inspect the active pPXF installation without modifying it.")
    mode.add_argument("--apply", action="store_true", help="Apply the pPXF 9.4.8 cached-whitener patch.")
    mode.add_argument("--restore", action="store_true", help="Restore ppxf.py from the backup made by --apply.")
    parser.add_argument(
        "--allow-outside-prefix",
        action="store_true",
        help="Allow patching a ppxf.py that is outside sys.prefix (not recommended for CRD_DAP).",
    )
    args = parser.parse_args()

    ppxf_pkg, ppxf_module = _import_ppxf()
    version = _ppxf_version(ppxf_pkg)
    path = _source_path(ppxf_module)
    source = path.read_text(encoding="utf-8")
    backup = Path(str(path) + BACKUP_SUFFIX)

    _status(path, source, version)

    if version != SUPPORTED_PPXF_VERSION:
        raise RuntimeError(
            f"This patch is intentionally locked to pPXF {SUPPORTED_PPXF_VERSION}; "
            f"the active environment has {version}."
        )
    if not args.allow_outside_prefix and not _is_inside(path, Path(sys.prefix)):
        raise RuntimeError(
            "The imported ppxf.py is outside sys.prefix. This can happen with PYTHONPATH, "
            "--user installs, or editable/shared installs. Refusing to modify it because the "
            "change may affect environments other than CRD_DAP."
        )

    if args.check:
        if PATCH_VERSION not in source:
            # A dry-run source-layout validation is useful before the user edits anything.
            _patched_source(source)
            print("Patchability check : PASS")
        else:
            print("Patchability check : already patched")
        return 0

    if args.restore:
        if not backup.is_file():
            raise FileNotFoundError(f"No CRD_DAP backup exists: {backup}")
        shutil.copy2(backup, path)
        importlib.invalidate_caches()
        print(f"Restored original ppxf.py from: {backup}")
        return 0

    # --apply
    if PATCH_VERSION in source:
        print("No action needed: CRD_DAP cached-whitener patch is already installed.")
        return 0
    if backup.exists():
        raise RuntimeError(
            f"Backup already exists: {backup}\n"
            "Refusing to overwrite it. Restore first or inspect the existing backup manually."
        )

    patched = _patched_source(source)
    shutil.copy2(path, backup)

    # Atomic write. If the subprocess import check fails, immediately restore.
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=str(path.parent), prefix=path.name + ".", suffix=".tmp", delete=False
        ) as handle:
            handle.write(patched)
            temp_path = Path(handle.name)
        os.replace(temp_path, path)
        importlib.invalidate_caches()
        _subprocess_smoke()
    except Exception:
        shutil.copy2(backup, path)
        importlib.invalidate_caches()
        raise

    print("Patch applied successfully.")
    print(f"Modified file : {path}")
    print(f"Backup        : {backup}")
    print("Next step     : run the cached-whitener pytest file in the CRD_DAP environment.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
