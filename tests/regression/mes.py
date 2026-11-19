"""Operator errors of the method of exact solutions (MES) cases.

Every case under ``MES/`` saves its error field ``e``, the numerical operator
minus the analytic one. This script runs every input of every case at two small
resolutions and compares the largest absolute error to a committed reference,
so a change of any operator or boundary condition shows up.

The cases save the numerical operator ``S_num`` and the error ``e = S_num - S``
after one evaluation of ``rhs()``, but the operators leave the guard cells
unset (zero), so there ``e = -S``. The largest error over all points was that
source in the outermost guard cell: independent of the operator, the method and
the resolution. The error is therefore taken over the inner points only, except
for the boundary condition cases under ``MES/boundaries``, which set the guard
cells and are measured there. The
script then proves it cannot go blind again: it reruns three cases with a
changed derivative method and requires the error to move far beyond the
tolerance.

Run from the repository root as ``python3 tests/regression/mes.py``. The
numbers of the run are written to ``recorded-mes.json`` in the current directory.
"""

import json
import math
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
# The volume integral of pure sines has the exact value zero, so its "error" is
# summation round-off of about 4e-12 that moves with the order of the sum
ATOL = 1.0e-12
RESOLUTIONS = ((32, 34, 64), (64, 66, 128))
# Case, changed options: each must move the error of the case
MUTATIONS = {
    "singleOperators/4-DDZCylinder/gaussianWSinAndParabola": ("mesh:ddz:first=C2",),
    "singleOperators/5-D2DZ2Cylinder/gaussianWSinAndParabola": ("mesh:ddz:second=C2",),
    "singleOperators/1-DDXCylinder/gaussian": ("mesh:ddx:first=W2",),
}
MUTATION_MARGIN = 1000.0


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


def inner(values: np.ndarray, mxg: int, myg: int) -> np.ndarray:
    """Cut the guard cells in x and y from a saved field.

    Parameters
    ----------
    values : np.ndarray
        The saved field, ``(x, y, z)`` or ``(t, x, y, z)`` (last step used);
        anything else (a scalar) is returned as it is.
    mxg : int
        Guard cells in x.
    myg : int
        Guard cells in y.

    Returns
    -------
    np.ndarray
        The inner points.
    """
    data = np.array(values)
    if data.ndim == 4:
        data = data[-1]
    if data.ndim != 3:
        return data
    return data[mxg : data.shape[0] - mxg, myg : data.shape[1] - myg, :]


def cases() -> List[Tuple[str, Path, Path]]:
    """Find every MES input.

    Returns
    -------
    list
        Tuples of the case name, the executable and the input file.
    """
    found = []
    for makefile in sorted((ROOT / "MES").rglob("CMakeLists.txt")):
        source = re.search(r"add_executable\((\w+)\s", makefile.read_text())
        if source is None:
            continue
        folder = makefile.parent
        for inp in sorted(folder.glob("*/BOUT.inp")):
            name = str(inp.parent.relative_to(ROOT / "MES"))
            found.append((name, folder / source.group(1), inp))
    return found


def run_case(
    exe: Path,
    inp: Path,
    resolution: Tuple[int, int, int],
    extra: Tuple[str, ...] = (),
) -> Optional[float]:
    """Run one input at one resolution.

    Parameters
    ----------
    exe : Path
        The compiled case.
    inp : Path
        Its input file.
    resolution : tuple
        Inner points in x, points in y and in z, each capped at the input's own.
    extra : tuple
        Further command line options, for the mutations.

    Returns
    -------
    float or None
        The largest absolute error, ``None`` if the run failed.
    """
    text = inp.read_text()
    mxg = read_int(text, "MXG", 1)
    myg = read_int(text, "MYG", 1)
    nxi, nyi, nzi = resolution
    own_ny = read_int(text, "ny", 1)
    own_nz = read_int(text, "nz", 1)
    nx = min(nxi, read_int(text, "nx", nxi + 2 * mxg) - 2 * mxg) + 2 * mxg
    ny = min(nyi, own_ny)
    nz = min(nzi, own_nz)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(inp, Path(tmp) / "BOUT.inp")
        command = [str(exe), "-d", tmp, f"mesh:nx={nx}", f"mesh:ny={ny}"]
        command += [f"mesh:nz={nz}", "NXPE=1", "input:error_on_unused_options=false"]
        command += list(extra)
        done = subprocess.run(  # nosec
            command,
            cwd=exe.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=600,
        )
        dump = Path(tmp) / "BOUT.dmp.0.nc"
        if done.returncode != 0 or not dump.exists():
            lines = [x for x in done.stdout.splitlines() if x.strip()]
            keys = ("Error encountered", "terminate called", "Assertion")
            first = next(
                (i for i, x in enumerate(lines) if any(k in x for k in keys)), -1
            )
            tail = lines[first : first + 8] if first >= 0 else lines[-6:]
            print(
                f"  run failed ({exe.name}, {resolution}): "
                + " | ".join(t[:300] for t in tail)
            )
            return None
        with Dataset(dump) as data:
            error = np.array(data.variables["e"][:])
            if "boundaries" not in inp.parts:
                error = inner(error, mxg, myg)
            return float(np.max(np.abs(error)))


def order(errors: List[Optional[float]]) -> str:
    """Observed order of convergence between the two resolutions.

    Parameters
    ----------
    errors : list
        Errors at the two resolutions.

    Returns
    -------
    str
        The order, or ``n/a``.
    """
    old, new = errors
    if not old or not new:
        return "n/a"
    return f"{math.log2(old / new):.2f}"


def check_mutations(reference: Dict[str, List[Optional[float]]]) -> bool:
    """Prove that the cases evaluate their operator, by changing it.

    Parameters
    ----------
    reference : dict
        The committed errors.

    Returns
    -------
    bool
        True if every mutation moved its case far beyond the tolerance.
    """
    paths = {name: (exe, inp) for name, exe, inp in cases()}
    good = True
    for name, extra in MUTATIONS.items():
        exe, inp = paths[name]
        mutated = run_case(exe, inp, RESOLUTIONS[0], extra)
        old = reference[name][0]
        if mutated is None or old is None:
            print(f"FAIL mutation {extra} of {name} did not run")
            good = False
            continue
        move = abs(mutated - old) / (RTOL * abs(old) + ATOL)
        status = "ok  " if move >= MUTATION_MARGIN else "FAIL"
        good = good and move >= MUTATION_MARGIN
        print(
            f"{status} mutation {extra} of {name}: error {old:.3e} -> "
            f"{mutated:.3e}, {move:.1e} x tolerance"
        )
    return good


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
        print(name, actual[name], "order", order(actual[name]))
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
                if not abs(new - old) <= RTOL * abs(old) + ATOL:
                    failed.append(name)
    for name in sorted(set(failed)):
        print(f"FAIL {name}: {actual.get(name)} against {reference.get(name)}")
    blind = not check_mutations(reference)
    return 1 if failed or blind else 0


if __name__ == "__main__":
    sys.exit(main())
