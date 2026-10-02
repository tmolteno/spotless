#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#
"""Stop-gap for GitHub issue tmolteno/spotless#1: --fits is not implemented yet.

Writing a FITS image needs a coordinate-aware sphere so that an image
array can be built to hand to tart_tools.api_imaging.save_fits_image();
that work is tracked in tmolteno/disko#10. Until it lands, every
handle_image() call site passes img=None, so requesting --fits must
report "not implemented" cleanly (message on stderr, non-zero exit)
instead of dying with an AttributeError traceback inside tart_tools.
"""

import sys

FITS_NOT_IMPLEMENTED_MSG = (
    "FITS export is not implemented yet (tracked in tmolteno/disko#10); "
    "--fits cannot write a file from this path"
)


def abort_fits_not_implemented():
    """Report the unimplemented --fits path on stderr and exit with status 2.

    Status 2 matches argparse's convention for user argument errors.
    No traceback is printed: this is a request for a feature that does
    not exist yet, not a crash.
    """
    print(f"error: {FITS_NOT_IMPLEMENTED_MSG}", file=sys.stderr)
    raise SystemExit(2)
