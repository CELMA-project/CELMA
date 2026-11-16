"""Operator errors of the method of exact solutions (MES) cases.

Every case under ``MES/`` saves its error field ``e``, the numerical operator
minus the analytic one. This script runs every input of every case at two small
resolutions and compares the largest absolute error to a committed reference,
so a change of any operator or boundary condition shows up.

Run from the repository root as ``python3 tests/regression/mes.py``. The
numbers of the run are written to ``recorded-mes.json`` in the current directory.
"""

import json
import re
import shutil
import subprocess  # nosec
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path(__file__).with_name("mes-reference.json")
RECORDED = Path("recorded-mes.json")
RTOL = 1.0e-6
NOISE = 1.0e-12
RESOLUTIONS = ((32, 34, 64), (64, 66, 128))


def read_int(text: str, key: str, default: int) -> int:
    """Read an integer option, also written as ``2^8``, from an input file.

    Parameters
    ----------
    text : str
        Content of the input file.
    key : str
        Option name.
    default : int
        Value if the option is not found.

    Returns
    -------
    int
        The value.
    """
    match = re.search(rf"^\s*{key}\s*=\s*(\d+)(?:\s*\^\s*(\d+))?", text, re.M)
    if not match:
        return default
    base, power = int(match.group(1)), match.group(2)
    return base ** int(power) if power else base


def cases() -> List[Tuple[str, Path, Path]]:
    """Find every MES input.

    Returns
    -------
    list
        Tuples of the case name, the executable and the input file.
    """
    found = []
    for makefile in sorted((ROOT / "MES").rglob("makefile")):
        source = re.search(r"SOURCEC\s*=\s*(\w+)\.cxx", makefile.read_text())
        if source is None:
            continue
        folder = makefile.parent
        for inp in sorted(folder.glob("*/BOUT.inp")):
            name = str(inp.parent.relative_to(ROOT / "MES"))
            found.append((name, folder / source.group(1), inp))
    return found


def run_case(exe: Path, inp: Path, resolution: Tuple[int, int, int]) -> Optional[float]:
    """Run one input at one resolution.

    Parameters
    ----------
    exe : Path
        The compiled case.
    inp : Path
        Its input file.
    resolution : tuple
        Inner points in x, points in y and in z, each capped at the input's own.

    Returns
    -------
    float or None
        The largest absolute error, ``None`` if the run failed.
    """
    text = inp.read_text()
    mxg = read_int(text, "MXG", 1)
    nxi, nyi, nzi = resolution
    own_ny = read_int(text, "ny", 1)
    own_nz = read_int(text, "nz", 1)
    nx = min(nxi, read_int(text, "nx", nxi + 2 * mxg) - 2 * mxg) + 2 * mxg
    ny = min(nyi, own_ny)
    nz = min(nzi, own_nz)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(inp, Path(tmp) / "BOUT.inp")
        command = [str(exe), "-d", tmp, f"mesh:nx={nx}", f"mesh:ny={ny}"]
        command += [f"mesh:nz={nz}", "NXPE=1"]
        done = subprocess.run(  # nosec
            command,
            cwd=exe.parent,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=600,
        )
        dump = Path(tmp) / "BOUT.dmp.0.nc"
        if done.returncode != 0 or not dump.exists():
            return None
        with Dataset(dump) as data:
            return float(np.max(np.abs(np.array(data.variables["e"][:]))))


def main() -> int:
    """Run every case and compare with the reference.

    Returns
    -------
    int
        Zero when every error agrees with the reference.
    """
    actual: Dict[str, List[Optional[float]]] = {}
    for name, exe, inp in cases():
        actual[name] = [run_case(exe, inp, res) for res in RESOLUTIONS]
        print(name, actual[name])
    RECORDED.write_text(json.dumps(actual, indent=1) + "\n")
    if not REFERENCE.exists():
        print(f"FAIL no reference at {REFERENCE}; the numbers are in {RECORDED}")
        return 1
    reference = json.loads(REFERENCE.read_text())
    failed = sorted(set(actual) ^ set(reference))
    for name in sorted(set(actual) & set(reference)):
        for new, old in zip(actual[name], reference[name]):
            if (new is None) != (old is None):
                failed.append(name)
            elif new is not None and old is not None:
                if not abs(new - old) <= RTOL * max(abs(old), NOISE):
                    failed.append(name)
    for name in sorted(set(failed)):
        print(f"FAIL {name}: {actual.get(name)} against {reference.get(name)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
