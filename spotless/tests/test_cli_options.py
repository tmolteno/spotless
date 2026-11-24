#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#
"""Options spotless no longer has, because disko_draw does the job now.

``--show-sources`` and ``--SVG`` were removed from ``spotless``: sources and SVG
are drawn by disko_draw, on the field of view that ``spotless --HDF`` writes. A
script that still passes one should be told where it went, not get argparse's
bare "unrecognized arguments", and ``--help`` should not advertise it or the
options that only existed to feed it.
"""

import unittest

from spotless.tests.test_fits_export import run_cli_script

CLIS = (("spotless", "spotless.spotless_cli"),
        ("gridless", "spotless.gridless_cli"))
# (option, the commands that used to take it)
REMOVED = (("--show-sources", CLIS), ("--SVG", CLIS[:1]))
GONE_FROM_HELP = ("--show-sources", "--SVG", "--elevation", "--catalog")


class TestRemovedOptions(unittest.TestCase):

    def test_the_old_options_point_at_disko_draw(self):
        for option, clis in REMOVED:
            for script, module in clis:
                with self.subTest(script=script, option=option):
                    proc = run_cli_script(script, module, [option])
                    self.assertNotEqual(
                        proc.returncode, 0, f"expected non-zero exit, stderr:\n{proc.stderr}"
                    )
                    self.assertIn(f"{option} has moved to disko_draw", proc.stderr)
                    self.assertNotIn("Traceback", proc.stderr)

    def test_the_options_are_no_longer_advertised(self):
        for script, module in CLIS:
            with self.subTest(script):
                proc = run_cli_script(script, module, ["--help"])
                self.assertEqual(proc.returncode, 0, proc.stderr)
                for option in GONE_FROM_HELP:
                    self.assertNotIn(option, proc.stdout)


if __name__ == "__main__":
    unittest.main()
