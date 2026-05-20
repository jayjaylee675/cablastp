"""Database configuration (port of dbconf.go).

The textual format on disk is `key:value` per line. This module reads/writes
that format and supports merging command-line overrides over a file-loaded
configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, asdict, replace
from typing import IO, Mapping, Optional, Set


@dataclass
class DBConf:
    MinMatchLen: int = 40
    MatchKmerSize: int = 4
    GappedWindowSize: int = 25
    UngappedWindowSize: int = 10
    ExtSeqIdThreshold: int = 60
    MatchSeqIdThreshold: int = 70
    MatchExtend: int = 30
    MapSeedSize: int = 6
    ExtSeedSize: int = 0
    LowComplexity: int = 10
    SeedLowComplexity: int = 6
    SavePlain: bool = False
    ReadOnly: bool = True
    BlastMakeBlastDB: str = "makeblastdb"
    BlastDBSize: int = 0

    def write(self, stream: IO) -> None:
        records = [
            ("MinMatchLen", str(self.MinMatchLen)),
            ("MatchKmerSize", str(self.MatchKmerSize)),
            ("GappedWindowSize", str(self.GappedWindowSize)),
            ("UngappedWindowSize", str(self.UngappedWindowSize)),
            ("ExtSeqIdThreshold", str(self.ExtSeqIdThreshold)),
            ("MatchSeqIdThreshold", str(self.MatchSeqIdThreshold)),
            ("MatchExtend", str(self.MatchExtend)),
            ("MapSeedSize", str(self.MapSeedSize)),
            ("ExtSeedSize", str(self.ExtSeedSize)),
            ("LowComplexity", str(self.LowComplexity)),
            ("SeedLowComplexity", str(self.SeedLowComplexity)),
            ("SavePlain", "1" if self.SavePlain else "0"),
            ("ReadOnly", "1" if self.ReadOnly else "0"),
            ("BlastMakeBlastDB", self.BlastMakeBlastDB),
            ("BlastDBSize", str(self.BlastDBSize)),
        ]
        out = "".join("%s:%s\n" % (k, v) for k, v in records)
        # Accept text- or binary-mode files.
        try:
            stream.write(out)
        except TypeError:
            stream.write(out.encode("utf-8"))

    Write = write

    def flag_merge(
        self,
        file_conf: "DBConf",
        explicit: Set[str],
    ) -> "DBConf":
        """Return a configuration where flags not in `explicit` are taken from
        `file_conf`. `explicit` is the set of flag names present on the command
        line (e.g., 'min-match-len').
        """
        if "map-seed-size" in explicit:
            raise ValueError(
                "The map seed size cannot be changed for an existing database."
            )
        if "read-only" in explicit:
            raise ValueError(
                "The read-only setting cannot be changed for an existing database."
            )

        merged = DBConf(**asdict(self))
        if "min-match-len" not in explicit:
            merged.MinMatchLen = file_conf.MinMatchLen
        if "match-kmer-size" not in explicit:
            merged.MatchKmerSize = file_conf.MatchKmerSize
        if "gapped-window-size" not in explicit:
            merged.GappedWindowSize = file_conf.GappedWindowSize
        if "ungapped-window-size" not in explicit:
            merged.UngappedWindowSize = file_conf.UngappedWindowSize
        if "ext-seq-id-threshold" not in explicit:
            merged.ExtSeqIdThreshold = file_conf.ExtSeqIdThreshold
        if "match-seq-id-threshold" not in explicit:
            merged.MatchSeqIdThreshold = file_conf.MatchSeqIdThreshold
        if "match-extend" not in explicit:
            merged.MatchExtend = file_conf.MatchExtend
        if "ext-seed-size" not in explicit:
            merged.ExtSeedSize = file_conf.ExtSeedSize
        if "low-complexity" not in explicit:
            merged.LowComplexity = file_conf.LowComplexity
        if "seed-low-complexity" not in explicit:
            merged.SeedLowComplexity = file_conf.SeedLowComplexity
        if "plain" not in explicit:
            merged.SavePlain = file_conf.SavePlain
        if "read-only" not in explicit:
            merged.ReadOnly = file_conf.ReadOnly
        if "makeblastdb" not in explicit:
            merged.BlastMakeBlastDB = file_conf.BlastMakeBlastDB
        if "dbsize" not in explicit:
            merged.BlastDBSize = file_conf.BlastDBSize
        return merged

    FlagMerge = flag_merge


# A snapshot of the default configuration. The Go code exposed a *pointer* and
# mutated it; in Python we hand out fresh copies via `default_db_conf()`.
DEFAULT_DB_CONF = DBConf()


def default_db_conf() -> DBConf:
    return DBConf(**asdict(DEFAULT_DB_CONF))


def load_db_conf(stream: IO) -> DBConf:
    """Load a configuration from a `key:value` text stream."""
    conf = default_db_conf()

    if hasattr(stream, "read"):
        data = stream.read()
        if isinstance(data, bytes):
            data = data.decode("utf-8")
    else:
        data = str(stream)

    for raw_line in data.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValueError("Invalid DBConf line: %r" % raw_line)
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if key == "MinMatchLen":
            conf.MinMatchLen = int(value)
        elif key == "MatchKmerSize":
            conf.MatchKmerSize = int(value)
        elif key == "GappedWindowSize":
            conf.GappedWindowSize = int(value)
        elif key == "UngappedWindowSize":
            conf.UngappedWindowSize = int(value)
        elif key == "ExtSeqIdThreshold":
            conf.ExtSeqIdThreshold = int(value)
        elif key == "MatchSeqIdThreshold":
            conf.MatchSeqIdThreshold = int(value)
        elif key == "MatchExtend":
            conf.MatchExtend = int(value)
        elif key == "MapSeedSize":
            conf.MapSeedSize = int(value)
        elif key == "ExtSeedSize":
            conf.ExtSeedSize = int(value)
        elif key == "LowComplexity":
            conf.LowComplexity = int(value)
        elif key == "SeedLowComplexity":
            conf.SeedLowComplexity = int(value)
        elif key == "SavePlain":
            conf.SavePlain = value == "1"
        elif key == "ReadOnly":
            conf.ReadOnly = value == "1"
        elif key == "BlastMakeBlastDB":
            conf.BlastMakeBlastDB = value
        elif key == "BlastDBSize":
            conf.BlastDBSize = int(value)
        else:
            raise ValueError("Invalid DBConf flag: %s" % key)

    return conf


LoadDBConf = load_db_conf
