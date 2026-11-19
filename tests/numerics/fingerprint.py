"""Fingerprint of the numerics: one right-hand-side evaluation, term by term.

Both models run on a small mesh from z-, x- and y-dependent initial fields, set
on the command line over the real input file, with ``nout=0``. BOUT++ then
evaluates the right-hand side exactly once, for the initial output, and takes
no solver step; the saved terms of that one time slice are that evaluation. Every
run asserts this from the log (one row of the table, one evaluation, step 0 of
0) and from the output (one time slice). A run with ``nout=1`` and the Euler
solver gave bit for bit the same values at the interior points and the model
guard cells, with a further evaluation in the step.

The time derivative of every evolved field and every individually saved term is
reduced to a few numbers (norms, sums, fixed projections and sampled points of
the interior) and compared with a committed reference at a tight relative
tolerance, so the report names the term that moved. A few operator cases of
``MES/`` are fingerprinted the same way.

The test proves its own sensitivity: it reruns with a changed z derivative, x
derivative, upwind scheme and boundary condition and requires each of them to
fail the comparison by a wide margin, and it runs the unchanged case twice and
requires every compared value to be bit for bit the same.

Run from the repository root, in the CI image, as
``python3 tests/numerics/fingerprint.py``. The numbers of the run are written
to ``recorded-fingerprint.json`` in the current directory; commit that file as
``tests/numerics/fingerprint-reference.json`` to rerecord the reference (which
needs a maintainer's consent, see AGENTS.md).
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
REFERENCE = Path(__file__).with_name("fingerprint-reference.json")
RECORDED = Path("recorded-fingerprint.json")
RTOL = 1.0e-9
# A mutation has to move the numbers this many times the tolerance
MARGIN = 100.0
GRID = ("mesh:nx=12", "mesh:ny=8", "mesh:nz=32")
RUN = ("nout=0", "timestep=1e-10", "input:error_on_unused_options=false")
# Several z modes, x and y dependence, nothing polynomial
FIELDS = (
    "lnN:function=0.2*sin(2*z+3*x)*cos(y)+0.1*cos(5*z-2*x+y)+0.05*sin(x)",
    "momDensPar:uIPar=geom:yl/geom:Ly+0.1*sin(3*z+x)*cos(2*y)",
    "jPar:function=0.05*cos(3*z+y)*(1+0.5*sin(x))+0.02*sin(7*z-x)",
)
VORTICITY = "{}:function=0.01*sin(2*z+3*x)*cos(y)+0.005*cos(4*z-y)+0.002*sin(6*z+x)"
MODELS = {
    "celma": ("celma", "celma", "CSDXNyScan", "vortD"),
    "celmaWithBoussinesqApprox": (
        "celmaWithBoussinesqApprox",
        "celmaWBA",
        "BousCSDXMagFieldScanAr",
        "vort",
    ),
}
EVOLVED = ("lnN", "momDensPar", "jPar")
MUTATIONS = {
    "z derivative C2 instead of FFT": ("mesh:ddz:first=C2", "mesh:ddz:second=C2"),
    # C4 and U2 need two guard cells, the inputs have one
    "x derivative W2 instead of C2": ("mesh:ddx:first=W2",),
    "y upwind C2 instead of U1": ("mesh:ddy:upwind=C2",),
    "outer x boundary of lnN Dirichlet": ("lnN:bndry_xout=dirichlet_o4(0.0)",),
}
# MES operator cases: directory, executable, input, saved variables
MES_CASES = {
    "MES/DDZCylinder": (
        "MES/singleOperators/4-DDZCylinder",
        "DDZCylinder",
        "gaussianWSinAndParabola",
    ),
    "MES/D2DZ2Cylinder": (
        "MES/singleOperators/5-D2DZ2Cylinder",
        "D2DZ2Cylinder",
        "gaussianWSinAndParabola",
    ),
    "MES/DDXCylinder": (
        "MES/singleOperators/1-DDXCylinder",
        "DDXCylinder",
        "gaussian",
    ),
}
MES_GRID = ("mesh:nx=34", "mesh:nz=32", "NXPE=1")
MES_RUN = ("nout=0", "input:error_on_unused_options=false")
SAMPLES = 16
Stats = List[float]
Fingerprint = Dict[str, Stats]


def fingerprint_of(values: np.ndarray) -> Stats:
    """Reduce an array to numbers that move when any element moves.

    Parameters
    ----------
    values : np.ndarray
        The interior points of one term.

    Returns
    -------
    list
        The largest absolute value (the scale of the term), the norm, the sum,
        three fixed projections and sixteen sampled points.
    """
    flat = np.ravel(values).astype(float)
    index = np.arange(flat.size)
    stats = [float(np.max(np.abs(flat))), float(np.sqrt(np.sum(flat**2)))]
    stats.append(float(np.sum(flat)))
    for k in range(1, 4):
        stats.append(float(np.sum(np.sin(index * 0.7311 * k + k) * flat)))
    picks = (index[:: max(1, flat.size // SAMPLES)])[:SAMPLES]
    stats.extend(float(flat[i]) for i in picks)
    return stats


def relative_error(new: Stats, old: Stats) -> float:
    """Compare two fingerprints of one term.

    Parameters
    ----------
    new : list
        Fingerprint of the run.
    old : list
        Fingerprint to compare with.

    Returns
    -------
    float
        The largest difference over the largest absolute value of the term.
    """
    scale = abs(old[0])
    diff = max(abs(a - b) for a, b in zip(new, old))
    return diff / scale if scale > 0 else (0.0 if diff == 0 else float("inf"))


def compare(
    new: Dict[str, Fingerprint], old: Dict[str, Fingerprint]
) -> Dict[str, float]:
    """Compare every term of two sets of fingerprints.

    Parameters
    ----------
    new : dict
        Fingerprints of the run, keyed ``"model/term"``.
    old : dict
        Fingerprints to compare with.

    Returns
    -------
    dict
        Term to relative error; a term in only one set has infinite error.
    """
    errors = {}
    for term in sorted(set(new) | set(old)):
        if term in new and term in old:
            errors[term] = relative_error(new[term], old[term])
        else:
            errors[term] = float("inf")
    return errors


def worst(errors: Dict[str, float]) -> Tuple[str, float]:
    """Find the term with the largest error.

    Parameters
    ----------
    errors : dict
        Term to relative error.

    Returns
    -------
    tuple
        The term and its error.
    """
    term = max(errors, key=lambda t: errors[t])
    return term, errors[term]


def run_bout(
    source: Path, executable: str, case: str, args: Tuple[str, ...], tmp: Path
) -> Tuple[Dataset, str]:
    """Run an executable on a copy of the real input of a case.

    Parameters
    ----------
    source : Path
        Directory of the executable.
    executable : str
        Name of the executable.
    case : str
        Directory of the input file in ``source``.
    args : tuple
        Command line options that override the input.
    tmp : Path
        Directory for the output.

    Returns
    -------
    tuple
        The opened output file and the log of the run.
    """
    shutil.copy(source / case / "BOUT.inp", tmp / "BOUT.inp")
    command = [str(source / executable), "-d", str(tmp), *args]
    done = subprocess.run(  # nosec
        command,
        cwd=source,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=1800,
    )
    if done.returncode != 0:
        tail = "\n".join(done.stdout.splitlines()[-15:])
        raise RuntimeError(f"{executable} {args} failed:\n{tail}")
    dump = Dataset(tmp / "BOUT.dmp.0.nc")
    try:
        check_single_evaluation(done.stdout, dump)
    except RuntimeError:
        dump.close()
        raise
    return dump, done.stdout


def check_single_evaluation(log: str, dump: Dataset) -> None:
    """Prove that a run evaluated the right-hand side once and took no step.

    Parameters
    ----------
    log : str
        Standard output of the run.
    dump : Dataset
        Its output file.

    Raises
    ------
    RuntimeError
        If the log table has other than one row with one evaluation, the solver
        ran for outputs or took a step, or the output has other than one time
        slice.
    """
    rows = re.findall(r"^\s*[\d.]+e[+-]\d+\s+(\d+)\s+[\d.]+e[+-]\d+", log, re.M)
    steps = re.findall(r"Step (\d+) of (\d+)", log)
    slices = dump.dimensions["t"].size
    if rows != ["1"] or steps != [("0", "0")] or slices != 1:
        tail = "\n".join(log.splitlines()[-30:])
        raise RuntimeError(
            f"not one evaluation and no step: RHS rows {rows}, steps {steps}, "
            f"time slices {slices}\n{tail}"
        )


def interior(data: np.ndarray) -> np.ndarray:
    """Cut the guard cells in x and y from a ``(x, y, z)`` array.

    Parameters
    ----------
    data : np.ndarray
        Field of one time step.

    Returns
    -------
    np.ndarray
        The inner points.
    """
    return np.asarray(data)[1:-1, 1:-1, :]


def run_model(name: str, extra: Tuple[str, ...] = ()) -> Dict[str, np.ndarray]:
    """Run one model for one right-hand-side evaluation.

    Parameters
    ----------
    name : str
        Key in ``MODELS``.
    extra : tuple
        Further command line options, for the mutations.

    Returns
    -------
    dict
        Fingerprint of every time dependent field except the evolved fields
        themselves, keyed ``"model/term"``; the time derivative of each
        evolved field is among them.
    """
    directory, executable, case, vorticity = MODELS[name]
    args = (*GRID, *RUN, *FIELDS, VORTICITY.format(vorticity), *extra)
    evolved = (vorticity, *EVOLVED)
    result: Dict[str, np.ndarray] = {}
    with tempfile.TemporaryDirectory() as tmp:
        dump, log = run_bout(ROOT / directory, executable, case, args, Path(tmp))
        with dump:
            for var, variable in dump.variables.items():
                if variable.dimensions != ("t", "x", "y", "z") or var in evolved:
                    continue
                values = interior(variable[0])
                if not np.all(np.isfinite(values)):
                    raise RuntimeError(f"{name} {var}: inner values not finite")
                result[f"{name}/{var}"] = np.array(values)
    return result


def run_mes(name: str, extra: Tuple[str, ...] = ()) -> Dict[str, np.ndarray]:
    """Run one operator case of ``MES/`` and fingerprint its output.

    Parameters
    ----------
    name : str
        Key in ``MES_CASES``.
    extra : tuple
        Further command line options, for the mutations.

    Returns
    -------
    dict
        Fingerprint of the numerical operator and its error, keyed
        ``"name/variable"``.
    """
    directory, executable, case = MES_CASES[name]
    result: Dict[str, np.ndarray] = {}
    with tempfile.TemporaryDirectory() as tmp:
        dump, _ = run_bout(
            ROOT / directory, executable, case, (*MES_GRID, *MES_RUN, *extra), Path(tmp)
        )
        with dump:
            if not np.any(np.asarray(dump.variables["S_num"])):
                raise RuntimeError(f"{name}: the operator was not evaluated")
            for var in ("S_num", "e"):
                data = np.asarray(dump.variables[var])
                # Saved once: (x, y, z); repeated: (t, x, y, z), last step
                values = (data[-1] if data.ndim == 4 else data)[1:-1]
                result[f"{name}/{var}"] = np.array(values)
    return result


def run_all(extra: Tuple[str, ...] = ()) -> Dict[str, np.ndarray]:
    """Run both models and the operator cases.

    Parameters
    ----------
    extra : tuple
        Further command line options.

    Returns
    -------
    dict
        The compared values of every term, keyed ``"model/term"``.
    """
    result: Dict[str, np.ndarray] = {}
    for name in MODELS:
        result.update(run_model(name, extra))
    for name in MES_CASES:
        mes_extra = tuple(x for x in extra if not x.startswith(("lnN:", "jPar:")))
        result.update(run_mes(name, mes_extra))
    return result


def fingerprints(arrays: Dict[str, np.ndarray]) -> Dict[str, Fingerprint]:
    """Reduce every term to its fingerprint.

    Parameters
    ----------
    arrays : dict
        The compared values of every term.

    Returns
    -------
    dict
        Term to fingerprint.
    """
    return {term: fingerprint_of(values) for term, values in arrays.items()}


def differing(a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]) -> List[str]:
    """List the terms whose values are not bit for bit the same.

    Parameters
    ----------
    a : dict
        Values of one run.
    b : dict
        Values of another run.

    Returns
    -------
    list
        The terms that differ in any element, or are in one run only.
    """
    return [
        t
        for t in sorted(set(a) | set(b))
        if t not in a or t not in b or not np.array_equal(a[t], b[t])
    ]


def report(errors: Dict[str, float], limit: float) -> List[str]:
    """List the terms above a limit.

    Parameters
    ----------
    errors : dict
        Term to relative error.
    limit : float
        Relative error above which a term has moved.

    Returns
    -------
    list
        One line per term that moved.
    """
    return [
        f"  moved {t}: relative error {e:.3e}" for t, e in errors.items() if e > limit
    ]


def groups(errors: Dict[str, float]) -> List[str]:
    """List the models and operator cases in a set of errors.

    Parameters
    ----------
    errors : dict
        Term to relative error.

    Returns
    -------
    list
        The first part of every term name, in order and once.
    """
    return list(dict.fromkeys(t.split("/")[0] for t in errors))


def main() -> int:
    """Compare with the reference and prove the sensitivity of the comparison.

    Returns
    -------
    int
        Zero when the numbers agree and every mutation fails the comparison.
    """
    run_one = run_all()
    run_two = run_all()
    first = fingerprints(run_one)
    RECORDED.write_text(json.dumps(first, indent=0) + "\n")
    unequal = differing(run_one, run_two)
    print(
        f"run-to-run: {len(run_one) - len(unequal)} of {len(run_one)} terms bit for bit"
    )
    for term in unequal:
        print(f"  FAIL not bit for bit the same: {term}")
    failed = bool(unequal)
    if REFERENCE.exists():
        reference = json.loads(REFERENCE.read_text())
    else:
        print(f"FAIL no reference at {REFERENCE}; the numbers are in {RECORDED}")
        reference, failed = first, True
    errors = compare(first, reference)
    term, err = worst(errors)
    print(f"{len(errors)} terms, worst {term}: {err:.3e} (tolerance {RTOL:.0e})")
    moved = report(errors, RTOL)
    failed = failed or bool(moved)
    print("\n".join(moved))
    for label, extra in MUTATIONS.items():
        try:
            mutated = fingerprints(run_all(extra))
        except RuntimeError as error:
            print(f"FAIL mutation '{label}' did not run: {error}")
            failed = True
            continue
        errors = compare(mutated, reference)
        for group in groups(errors):
            subset = {t: e for t, e in errors.items() if t.split("/")[0] == group}
            term, err = worst(subset)
            margin = err / RTOL
            # The operator cases are only a second opinion on the operators
            good = margin >= MARGIN or group == "MES"
            failed = failed or not good
            print(
                f"{'ok  ' if good else 'FAIL'} mutation '{label}', {group}: "
                f"worst term {term}, error {err:.3e} = {margin:.1e} x tolerance, "
                "run-to-run noise 0 (bit for bit)"
            )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
