"""Throwaway probe: is rhs evaluated for nout=0 and nout=1?"""

import re
import shutil
import subprocess  # nosec
import tempfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

ROOT = Path(__file__).resolve().parents[2]
CASES = {
    "DDZ": ("MES/singleOperators/4-DDZCylinder", "DDZCylinder", "gaussianWSinAndParabola"),
    "D2DZ2": ("MES/singleOperators/5-D2DZ2Cylinder", "D2DZ2Cylinder", "gaussianWSinAndParabola"),
    "DDX": ("MES/singleOperators/1-DDXCylinder", "DDXCylinder", "gaussian"),
}
for name, (d, exe, case) in CASES.items():
    for nout in (0, 1):
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(ROOT / d / case / "BOUT.inp", Path(tmp) / "BOUT.inp")
            done = subprocess.run(  # nosec
                [str(ROOT / d / exe), "-d", tmp, "mesh:nx=34", "mesh:nz=32", "NXPE=1",
                 f"nout={nout}", "input:error_on_unused_options=false"],
                cwd=ROOT / d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
            rows = re.findall(r"^\s*[\d.]+e[+-]\d+\s+(\d+)\s+[\d.]+e[+-]\d+", done.stdout, re.M)
            with Dataset(Path(tmp) / "BOUT.dmp.0.nc") as dump:
                s_num, e, s = (np.asarray(dump.variables[v]) for v in ("S_num", "e", "S"))
            print(f"PROBE {name} nout={nout} rc={done.returncode} rhs_rows={rows} "
                  f"max|S_num|={np.max(np.abs(s_num)):.3e} max|e|={np.max(np.abs(e)):.3e} "
                  f"max|e+S|={np.max(np.abs(e + s)):.3e} max|S|={np.max(np.abs(s)):.3e}")
