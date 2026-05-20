"""DB: top-level cablastp database (port of db.go)."""

from __future__ import annotations

import os
import shutil
from typing import Optional, Set

from cablastp.coarse import CoarseDB, FILE_COARSE_FASTA
from cablastp.compressed import CompressedDB
from cablastp.dbconf import DBConf, default_db_conf, load_db_conf
from cablastp.exec_util import run_command
from cablastp.misc import vprintf

FILE_PARAMS = "params"
FILE_BLAST_COARSE = "blastdb-coarse"
FILE_BLAST_FINE = "blastdb-fine"


class DB:
    """A cablastp database directory: coarse + compressed parts plus config."""

    def __init__(self, path: str, conf: Optional[DBConf] = None):
        self.path = path
        self.name = os.path.basename(os.path.normpath(path))
        self.conf: DBConf = conf if conf is not None else default_db_conf()
        self.com_db: Optional[CompressedDB] = None
        self.coarse_db: Optional[CoarseDB] = None
        self._appending: bool = False
        self._params_file = None

    # ----------------------------------------------- DBConf delegation

    @property
    def MinMatchLen(self) -> int: return self.conf.MinMatchLen
    @property
    def MatchKmerSize(self) -> int: return self.conf.MatchKmerSize
    @property
    def GappedWindowSize(self) -> int: return self.conf.GappedWindowSize
    @property
    def UngappedWindowSize(self) -> int: return self.conf.UngappedWindowSize
    @property
    def ExtSeqIdThreshold(self) -> int: return self.conf.ExtSeqIdThreshold
    @property
    def MatchSeqIdThreshold(self) -> int: return self.conf.MatchSeqIdThreshold
    @property
    def MatchExtend(self) -> int: return self.conf.MatchExtend
    @property
    def MapSeedSize(self) -> int: return self.conf.MapSeedSize
    @property
    def ExtSeedSize(self) -> int: return self.conf.ExtSeedSize
    @property
    def LowComplexity(self) -> int: return self.conf.LowComplexity
    @property
    def SeedLowComplexity(self) -> int: return self.conf.SeedLowComplexity
    @property
    def SavePlain(self) -> bool: return self.conf.SavePlain
    @property
    def ReadOnly(self) -> bool: return self.conf.ReadOnly
    @property
    def BlastMakeBlastDB(self) -> str: return self.conf.BlastMakeBlastDB
    @property
    def BlastDBSize(self) -> int: return self.conf.BlastDBSize

    @BlastDBSize.setter
    def BlastDBSize(self, value: int) -> None:
        self.conf.BlastDBSize = value

    # Convenience aliases mirroring the embedded *DBConf in Go.
    @property
    def CoarseDB(self) -> CoarseDB:
        return self.coarse_db

    @property
    def ComDB(self) -> CompressedDB:
        return self.com_db

    @property
    def Path(self) -> str:
        return self.path

    @property
    def Name(self) -> str:
        return self.name

    # --------------------------------------------------- file utilities

    def file_path(self, name: str) -> str:
        return os.path.join(self.path, name)

    def open_write_file(self, append: bool, name: str):
        path = os.path.join(self.path, name)
        if append:
            return open(path, "r+b")
        return open(path, "w+b")

    def open_read_file(self, name: str):
        return open(os.path.join(self.path, name), "rb")

    # -------------------------------------------------------- save/close

    def save(self) -> None:
        if not self._appending:
            self._params_file.truncate(0)
            self._params_file.seek(0, os.SEEK_SET)
            self.conf.write(self._params_file)
            self._params_file.flush()

        self.coarse_db.save()

        # Build the coarse blastp database with makeblastdb.
        vprintf("Creating %s...\n", FILE_BLAST_COARSE)
        run_command(
            [self.BlastMakeBlastDB, "-dbtype", "prot",
             "-in", FILE_COARSE_FASTA, "-out", FILE_BLAST_COARSE],
            cwd=self.path,
        )
        vprintf("Done creating %s.\n", FILE_BLAST_COARSE)

    Save = save

    def read_close(self) -> None:
        if self._params_file:
            self._params_file.close()
        self.coarse_db.read_close()
        self.com_db.read_close()

    ReadClose = read_close

    def write_close(self) -> None:
        if self._params_file:
            self._params_file.close()
        self.coarse_db.write_close()
        self.com_db.write_close()

    WriteClose = write_close


def new_write_db(
    append: bool,
    conf: DBConf,
    directory: str,
    explicit_flags: Optional[Set[str]] = None,
) -> DB:
    """Create a database for writing, or open an existing one for appending.

    `explicit_flags` is the set of CLI flag names that were explicitly set by
    the user (e.g., {'min-match-len'}); used to merge with on-disk parameters
    when appending.
    """
    vprintf("Opening database in %s...\n", directory)

    if directory.endswith(".tar") or directory.endswith(".gz"):
        raise ValueError(
            "The CaBLASTP database you've provided does not appear to be a "
            "directory. Please make sure you've extracted the downloaded "
            "database with `tar zxf cablastp-xxx.tar.gz` before using it "
            "with CaBLASTP."
        )

    exists = os.path.isdir(directory)
    if append:
        if not exists:
            raise FileNotFoundError(
                "Could not open '%s' for appending: directory does not exist."
                % directory
            )
    else:
        if exists:
            raise FileExistsError(
                "The directory '%s' already exists. A new compressed database "
                "cannot be created in the same directory as an existing "
                "database. If you want to append to an existing database, "
                "use the '--append' flag." % directory
            )
        os.makedirs(directory, exist_ok=False)

    db = DB(directory, conf)
    db._appending = append

    # Sanity-check makeblastdb early so we don't do work for nothing.
    if shutil.which(db.BlastMakeBlastDB) is None:
        raise FileNotFoundError(
            "Could not find 'makeblastdb' executable: %s" % db.BlastMakeBlastDB
        )

    db._params_file = db.open_write_file(append, FILE_PARAMS)

    if append:
        db._params_file.seek(0, os.SEEK_SET)
        param_conf = load_db_conf(db._params_file)
        db.conf = db.conf.flag_merge(param_conf, explicit_flags or set())
        if db.conf.ReadOnly:
            raise ValueError("Appending to a read-only database is not possible.")

    db.com_db = CompressedDB.open_for_write(append, db)
    db.coarse_db = CoarseDB.open_for_write(append, db)

    vprintf("Done opening database in %s.\n", directory)
    return db


def new_read_db(directory: str) -> DB:
    """Open a cablastp database for reading."""
    vprintf("Opening database in %s...\n", directory)

    if directory.endswith(".tar") or directory.endswith(".gz"):
        raise ValueError(
            "The CaBLASTP database you've provided does not appear to be a "
            "directory. Please make sure you've extracted the downloaded "
            "database with `tar zxf cablastp-xxx.tar.gz` before using it "
            "with CaBLASTP."
        )

    if not os.path.isdir(directory):
        raise FileNotFoundError(
            "Could not open '%s' for reading: directory does not exist."
            % directory
        )

    db = DB(directory)
    db._params_file = db.open_read_file(FILE_PARAMS)
    db.conf = load_db_conf(db._params_file)

    if shutil.which(db.BlastMakeBlastDB) is None:
        raise FileNotFoundError(
            "Could not find 'makeblastdb' executable: %s" % db.BlastMakeBlastDB
        )

    db.com_db = CompressedDB.open_for_read(db)
    db.coarse_db = CoarseDB.open_for_read(db)

    vprintf("Done opening database in %s.\n", directory)
    return db


# Aliases in the Go casing.
NewWriteDB = new_write_db
NewReadDB = new_read_db
