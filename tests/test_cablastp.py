"""Pytest equivalents of cablastp_test.go."""

from __future__ import annotations

import io
from dataclasses import asdict

from cablastp.dbconf import DBConf, default_db_conf, load_db_conf
from cablastp.seeds import Seeds
from cablastp.seqdiff import (
    EditScript,
    new_edit_script,
    parse_edit_script,
)


def test_seed_hashing():
    seed_size = 4
    seed_low_complexity = 6
    seeds = Seeds(seed_size, seed_low_complexity)

    cases = [
        (b"ABCD", 578),
        (b"AAAA", 0),
    ]
    for kmer, expected_hash in cases:
        h = seeds.hash_kmer(kmer)
        assert h == expected_hash, (
            "kmer %r hashed to %d but should have hashed to %d"
            % (kmer, h, expected_hash)
        )
        recovered = seeds.unhash_kmer(expected_hash)
        assert recovered == kmer, (
            "hash %d unhashed to %r but should have unhashed to %r"
            % (expected_hash, recovered, kmer)
        )


def test_db_conf_io():
    conf = default_db_conf()
    buf = io.StringIO()
    conf.write(buf)

    buf.seek(0)
    conf_loaded = load_db_conf(buf)

    assert asdict(conf) == asdict(conf_loaded)


def test_edit_scripts():
    cases = [
        (
            # From the CaBLAST paper (nucleotide example, but the algorithm is
            # alphabet-agnostic so it still works as a sanity check).
            b"GTTCACTTATGTATTCATATGATTTTGGCAA",
            b"GTTCACGTGTATATTTATATAATTTTGGCAA",
            b"GTTCACTTATGTATTC--ATATGATTTTGGCAA",
            b"GTTCACG--TGTATATTTATATAATTTTGGCAA",
            "s6Gd1--s7ATi2TTs4A",
        ),
    ]
    for from_seq, to_seq, aligned_from, aligned_to, expected_script in cases:
        script = new_edit_script((aligned_from, aligned_to))
        assert str(script) == expected_script

        parsed = parse_edit_script(expected_script)
        assert str(parsed) == expected_script

        result = script.apply(from_seq)
        assert result == to_seq
