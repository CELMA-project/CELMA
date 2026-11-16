"""Regression test of the right-hand sides of the CELMA models.

Runs both models on a small mesh, collects the time derivatives and the
evolved fields from the output file and compares them to a committed reference.
The reference is recorded inside the continuous-integration image, as
the numbers depend on the BOUT++ version and the compiler.

Run from the repository root as ``python3 tests/regression/regression.py``.
"""

import gzip
import json
import shutil
import subprocess  # nosec
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path(__file__).with_name("reference.json.gz")
MODELS = {
    "celma": ("celma", "celma", "CSDXNyScan", "vortD"),
    "celmaWithBoussinesqApprox": (
        "celmaWithBoussinesqApprox",
        "celmaWBA",
        "BousCSDXMagFieldScanAr",
        "vort",
    ),
}
MESH = ("mesh:nx=12", "mesh:ny=8", "mesh:nz=32")
RUN = ("nout=10", "timestep=1e-2")
RTOL = 1.0e-6
ATOL = 1.0e-14
RECORDED = Path("recorded-reference.json")


def run_model(name: str, workdir: Path) -> Dict[str, List[float]]:
    """Run one model and collect the compared fields.

    Parameters
    ----------
    name : str
        Key in ``MODELS``.
    workdir : Path
        Directory in which the output is written.

    Returns
    -------
    dict
        Field name to the flattened values of the last time step.

    Raises
    ------
    KeyError
        If the output file lacks a compared field.
    """
    directory, executable, case, vorticity = MODELS[name]
    source = ROOT / directory
    shutil.copy(source / case / "BOUT.inp", workdir / "BOUT.inp")
    command = [str(source / executable), "-d", str(workdir), *MESH, *RUN]
    subprocess.run(command, cwd=source, check=True, timeout=1800)  # nosec
    wanted = ["ddt(" + v + ")" for v in (vorticity, "lnN", "momDensPar", "jPar")]
    wanted += [vorticity, "lnN", "momDensPar", "jPar"]
    fields = {}
    with Dataset(workdir / "BOUT.dmp.0.nc") as dump:
        missing = [w for w in wanted if w not in dump.variables]
        if missing:
            raise KeyError(f"{name} did not save {missing}")
        for field in wanted:
            fields[field] = np.ravel(dump.variables[field][-1]).tolist()
    return fields


def main() -> int:
    """Run both models and compare with the reference.

    Returns
    -------
    int
        Zero when every field agrees with the reference.
    """
    actual = {}
    for name in MODELS:
        with tempfile.TemporaryDirectory() as tmp:
            actual[name] = run_model(name, Path(tmp))
    for name, fields in actual.items():
        for field, values in fields.items():
            if not np.all(np.isfinite(values)):
                print(f"FAIL {name} {field}: not finite")
                return 1
    RECORDED.write_text(json.dumps(actual, indent=1) + "\n")
    if not REFERENCE.exists():
        print(f"FAIL no reference at {REFERENCE}; the numbers are in {RECORDED}")
        return 1
    reference = json.loads(gzip.decompress(REFERENCE.read_bytes()))
    failed = False
    for name, fields in actual.items():
        for field, values in fields.items():
            expected = np.array(reference[name][field])
            diff = float(np.max(np.abs(np.array(values) - expected)))
            allowed = RTOL * float(np.max(np.abs(expected))) + ATOL
            error = diff / allowed * RTOL
            status = "ok  " if error <= RTOL else "FAIL"
            failed = failed or not error <= RTOL
            print(f"{status} {name} {field}: relative error {error:.3e}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
