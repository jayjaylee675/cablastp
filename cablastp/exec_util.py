"""Subprocess wrapper that surfaces stderr in raised errors (port of cmd.go)."""

from __future__ import annotations

import subprocess
from typing import Optional, Sequence, Union, IO

from cablastp.misc import vprintf


def run_command(
    args: Sequence[str],
    stdin: Optional[Union[bytes, IO]] = None,
    stdout: Optional[IO] = None,
    cwd: Optional[str] = None,
) -> None:
    """Run a command and raise an error that includes stderr if it fails."""
    full_cmd = " ".join(args)
    vprintf("%s\n", full_cmd)

    stdin_arg = subprocess.PIPE if isinstance(stdin, (bytes, bytearray)) else stdin
    proc = subprocess.Popen(
        list(args),
        stdin=stdin_arg,
        stdout=stdout,
        stderr=subprocess.PIPE,
        cwd=cwd,
    )

    if isinstance(stdin, (bytes, bytearray)):
        _, stderr = proc.communicate(input=bytes(stdin))
    else:
        _, stderr = proc.communicate()

    if proc.returncode != 0:
        if stderr:
            raise RuntimeError(
                "Error running '%s': exit %d.\n\nstderr:\n%s"
                % (full_cmd, proc.returncode, stderr.decode("utf-8", errors="replace"))
            )
        raise RuntimeError(
            "Error running '%s': exit %d." % (full_cmd, proc.returncode)
        )
