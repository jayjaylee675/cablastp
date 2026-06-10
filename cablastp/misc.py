"""Verbose logging and flag formatting helpers (port of misc.go)."""

from __future__ import annotations

import sys
from typing import Any

Verbose: bool = False


def set_verbose(value: bool) -> None:
    global Verbose
    Verbose = value


def vprint(s: str) -> None:
    if not Verbose:
        return
    sys.stderr.write(s)


def vprintf(fmt: str, *args: Any) -> None:
    if not Verbose:
        return
    sys.stderr.write(fmt % args if args else fmt)


def vprintln(s: str) -> None:
    if not Verbose:
        return
    sys.stderr.write(s + "\n")


def print_flag_defaults(parser) -> None:
    """Print argparse flag defaults in the same style the Go CLI used."""
    for action in parser._actions:
        if not action.option_strings:
            continue
        # Pick the longest option string (the canonical form).
        name = max(action.option_strings, key=len).lstrip("-")
        default = action.default
        if default is None:
            default = ""
        elif isinstance(default, bool):
            # Go printed bool flag defaults lowercase ("true"/"false"), whereas
            # Python's str(True) is "True".
            default = str(default).lower()
        help_text = action.help or ""
        print('--%s="%s"\n\t%s' % (name, default, help_text))
