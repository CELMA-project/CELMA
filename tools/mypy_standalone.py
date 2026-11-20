"""Check independent scripts separately so duplicate module names cannot collide."""

import argparse
import subprocess
import sys


def main() -> int:
    """Run mypy for every supplied file and fail if any invocation fails."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ignore-missing-imports", action="store_true")
    parser.add_argument("--scripts-are-modules", action="store_true")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    options = []
    if args.ignore_missing_imports:
        options.append("--ignore-missing-imports")
    if args.scripts_are_modules:
        options.append("--scripts-are-modules")
    status = 0
    for filename in args.files:
        result = subprocess.run(
            [sys.executable, "-m", "mypy", *options, "--", filename], check=False
        )
        if result.returncode:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
