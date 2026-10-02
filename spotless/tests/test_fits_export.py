#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#
"""Regression tests for GitHub issue tmolteno/spotless#1.

--fits used to reach tart_tools.api_imaging.save_fits_image() with
img=None (every handle_image() call site passes None) and die with
``AttributeError: 'NoneType' object has no attribute 'astype'``. It must
now fail with a clear message naming the blocker (tmolteno/disko#10) on
stderr and a non-zero exit, with no Python traceback, for both CLIs.
"""

import io
import os
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from types import SimpleNamespace

from spotless.fits_export import FITS_NOT_IMPLEMENTED_MSG
from spotless.gridless_cli import handle_image as gridless_handle_image
from spotless.spotless_cli import handle_image as spotless_handle_image

# The venv's console scripts live next to the interpreter running pytest.
SCRIPT_DIR = os.path.dirname(sys.executable)


def run_cli_script(script, script_module, argv, timeout=300):
    """Run one of the CLI entry points as a real subprocess.

    Falls back to ``python -c`` if the console script is not installed
    (python puts '-c' in argv[0], so argparse still sees argv[1:]).
    """
    path = os.path.join(SCRIPT_DIR, script)
    if os.path.exists(path):
        cmd = [path]
    else:
        cmd = [
            sys.executable,
            "-c",
            f"from {script_module} import main; main()",
        ]
    return subprocess.run(
        cmd + argv,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


class TestFitsNotImplemented(unittest.TestCase):
    """--fits must fail cleanly (no traceback) while disko#10 is open."""

    def assert_clean_fits_failure(self, returncode, stderr):
        self.assertNotEqual(
            returncode, 0, f"expected non-zero exit, stderr:\n{stderr}"
        )
        self.assertIn(FITS_NOT_IMPLEMENTED_MSG, stderr)
        self.assertIn("tmolteno/disko#10", stderr)
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("AttributeError", stderr)

    def test_spotless_cli_fits_fails_cleanly(self):
        # The configuration from the bug report: --fits with data available.
        proc = run_cli_script(
            "spotless",
            "spotless.spotless_cli",
            [
                "--ms",
                "test_data/test.ms",
                "--fits",
                "--healpix",
                "--res",
                "120arcmin",
            ],
        )
        self.assert_clean_fits_failure(proc.returncode, proc.stderr)

    def test_gridless_cli_fits_fails_cleanly(self):
        proc = run_cli_script(
            "gridless",
            "spotless.gridless_cli",
            ["--file", "test_data/test_data.json", "--fits"],
        )
        self.assert_clean_fits_failure(proc.returncode, proc.stderr)

    def test_spotless_handle_image_without_image(self):
        # Directly exercise the crash site: --fits requested, img is None.
        args = SimpleNamespace(
            fits=True,
            PNG=False,
            PDF=False,
            SVG=False,
            display=False,
            dir=".",
            title="spotless",
        )
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as ctx:
                spotless_handle_image(args, None, "gridless", "2026_01_01_UTC")
        self.assertNotEqual(ctx.exception.code, 0)
        self.assertIn(FITS_NOT_IMPLEMENTED_MSG, stderr.getvalue())

    def test_gridless_handle_image_without_image(self):
        args = SimpleNamespace(
            fits=True,
            PNG=False,
            PDF=False,
            display=False,
            dir=".",
            title="",
        )
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as ctx:
                gridless_handle_image(args, None, "spotless", "2026_01_01_UTC")
        self.assertNotEqual(ctx.exception.code, 0)
        self.assertIn(FITS_NOT_IMPLEMENTED_MSG, stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
