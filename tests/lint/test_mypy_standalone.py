"""Exercise the standalone mypy command without running checked scripts."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

WRAPPER = Path(__file__).resolve().parents[2] / "tools" / "mypy_standalone.py"


class StandaloneMypyTests(unittest.TestCase):
    """Check isolation, diagnostics, and configuration at the CLI boundary."""

    def test_duplicate_names_and_spaces(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = [root / folder / "driver.py" for folder in ("first", "with space")]
            for path in paths:
                path.parent.mkdir()
                path.write_text('value: int = 1\nraise RuntimeError("never execute")\n')
            result = subprocess.run(
                [sys.executable, str(WRAPPER), *map(str, paths)],
                cwd=root,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_reports_every_failure_and_preserves_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "setup.cfg"
            config.write_text("[mypy]\ndisallow_untyped_calls=True\n")
            paths = [root / name for name in ("first.py", "last.py")]
            for path in paths:
                path.write_text("def untyped():\n    return 1\n\nuntyped()\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(WRAPPER),
                    "--ignore-missing-imports",
                    "--scripts-are-modules",
                    *map(str, paths),
                ],
                cwd=root,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            for path in paths:
                self.assertIn(path.name, result.stdout)
            self.assertIn("no-untyped-call", result.stdout)

    def test_no_files_is_an_error(self) -> None:
        result = subprocess.run(
            [sys.executable, str(WRAPPER)], capture_output=True, check=False
        )
        self.assertNotEqual(result.returncode, 0)

    def test_missing_import_option_is_forwarded(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "standalone.py"
            path.write_text("import celma_intentionally_missing_test_dependency\n")
            command = [sys.executable, str(WRAPPER), str(path)]
            missing = subprocess.run(
                command, cwd=directory, capture_output=True, check=False
            )
            ignored = subprocess.run(
                [*command[:2], "--ignore-missing-imports", *command[2:]],
                cwd=directory,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(missing.returncode, 0)
            self.assertEqual(ignored.returncode, 0, ignored.stdout + ignored.stderr)
