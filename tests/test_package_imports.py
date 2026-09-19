"""Checks for lightweight, workflow-specific package imports."""

import subprocess
import sys


def test_fuv_detector_does_not_import_unrelated_workflows():
    command = (
        "import sys; "
        "import icbuilder.fuvdetector; "
        "unexpected = {'secsy', 'icreader', 'icphysics', 'sksparse'} "
        "& set(sys.modules); "
        "assert not unexpected, sorted(unexpected)"
    )
    result = subprocess.run(
        [sys.executable, "-c", command],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
