"""Parse BLAST XML output (port of cmd/cablastp-*/xml.go).

The Go code used encoding/xml struct tags. Here we use ElementTree to read
out the same nested fields.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import IO, List, Union


@dataclass
class Hsp:
    num: int = 0
    evalue: float = 0.0
    query_from: int = 0
    query_to: int = 0
    hit_from: int = 0
    hit_to: int = 0


@dataclass
class Hit:
    num: int = 0
    accession: int = 0
    hsps: List[Hsp] = field(default_factory=list)


@dataclass
class BlastResults:
    hits: List[Hit] = field(default_factory=list)


def _text(elem, name: str, default: str = "") -> str:
    child = elem.find(name)
    if child is None or child.text is None:
        return default
    return child.text


def parse_blast_xml(source: Union[bytes, str, IO]) -> BlastResults:
    if isinstance(source, (bytes, bytearray)):
        root = ET.fromstring(bytes(source))
    elif isinstance(source, str):
        root = ET.fromstring(source)
    else:
        tree = ET.parse(source)
        root = tree.getroot()

    results = BlastResults()
    for hit_elem in root.iter("Hit"):
        hit = Hit(
            num=int(_text(hit_elem, "Hit_num", "0") or 0),
            accession=int(_text(hit_elem, "Hit_accession", "0") or 0),
        )
        hsps_container = hit_elem.find("Hit_hsps")
        if hsps_container is not None:
            for hsp_elem in hsps_container.findall("Hsp"):
                hsp = Hsp(
                    num=int(_text(hsp_elem, "Hsp_num", "0") or 0),
                    evalue=float(_text(hsp_elem, "Hsp_evalue", "0") or 0.0),
                    query_from=int(_text(hsp_elem, "Hsp_query-from", "0") or 0),
                    query_to=int(_text(hsp_elem, "Hsp_query-to", "0") or 0),
                    hit_from=int(_text(hsp_elem, "Hsp_hit-from", "0") or 0),
                    hit_to=int(_text(hsp_elem, "Hsp_hit-to", "0") or 0),
                )
                hit.hsps.append(hsp)
        results.hits.append(hit)
    return results
