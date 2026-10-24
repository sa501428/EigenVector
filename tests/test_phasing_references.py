#!/usr/bin/env python3
"""Independent sequence fixtures for reference building and rebinning correctness."""
import gzip
import hashlib
import json
import math
from argparse import Namespace
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from build_phasing_references import acquire, read_refseq
from export_phasing_track import gc_bins, refseq_bins
from phasing_reference import (TwoBit, count_chromosomes, gzip_output,
                               write_count_header)
from validate_phasing_references import independent_raw_gc, independently_decode


def fixture_2bit(path, bases, order="<"):
    def pack(*values):
        return struct.pack(order + "I" * len(values), *values)
    nblocks, masks = [], []
    for target, predicate in ((nblocks, lambda c: c.upper() == "N"), (masks, str.islower)):
        start = None
        for i, base in enumerate(bases + "!"):
            if predicate(base) and start is None:
                start = i
            if not predicate(base) and start is not None:
                target.append((start, i-start))
                start = None
    packed = bytearray()
    mapping = {"T": 0, "C": 1, "A": 2, "G": 3, "N": 3}  # encode Ns as G intentionally
    for start in range(0, len(bases), 4):
        byte = 0
        for i, base in enumerate(bases[start:start+4]):
            byte |= mapping[base.upper()] << (6-2*i)
        packed.append(byte)
    offset = 16 + 1 + 4 + 4
    data = pack(0x1A412743, 0, 1, 0) + b"\x04chr1" + pack(offset)
    data += pack(len(bases), len(nblocks))
    data += pack(*(s for s, n in nblocks)) + pack(*(n for s, n in nblocks))
    data += pack(len(masks)) + pack(*(s for s, n in masks)) + pack(*(n for s, n in masks)) + pack(0)
    path.write_bytes(data + packed)


def literal_counts(bases):
    bases = bases.upper()
    return bases.count("G") + bases.count("C"), sum(bases.count(c) for c in "ACGT")


class PhasingReferences(unittest.TestCase):
    def test_twobit_endianness_nblocks_and_partial_bytes(self):
        rng = random.Random(17)
        bases = "NNN" + "".join(rng.choice("ACGTacgt") for _ in range(512)) + "N"*111 + "gcATN"
        with tempfile.TemporaryDirectory() as directory:
            for endian in ("<", ">"):
                path = Path(directory) / "fixture.2bit"
                fixture_2bit(path, bases, endian)
                with TwoBit(path) as sequence:
                    for width in (1, 3, 4, 7, 25, 100, 125, 250, 1000):
                        records = list(sequence.counts("chr1", width))
                        self.assertEqual(sum(gc for _, _, gc, _ in records), literal_counts(bases)[0])
                        for start, end, gc, valid in records:
                            self.assertEqual((gc, valid), literal_counts(bases[start:end]))
                            self.assertEqual(independently_decode(sequence, "chr1", start, end), (gc, valid))
                    _, blocks, at = sequence.record("chr1")
                    total = independent_raw_gc(sequence, at, 0, len(bases))
                    total -= sum(independent_raw_gc(sequence, at, a, b) for a, b in blocks)
                    self.assertEqual(total, literal_counts(bases)[0])

    def test_count_roundtrip_and_terminal_tile(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.bin.gz"
            with gzip_output(path) as stream:
                write_count_header(stream, {"chr1": 201, "chr2": 1})
                stream.write(bytes([40, 100, 0, 0, 1, 1, 0, 1]))
            records = list(count_chromosomes(path))
            self.assertEqual(records, [("chr1", 201, bytes([40,100,0,0,1,1])), ("chr2", 1, bytes([0,1]))])
            first = path.read_bytes()
            with gzip_output(path) as stream:
                write_count_header(stream, {"chr1": 201, "chr2": 1})
                stream.write(bytes([40, 100, 0, 0, 1, 1, 0, 1]))
            self.assertEqual(first, path.read_bytes())

    def test_rebin_sums_base_counts_instead_of_percentages(self):
        data = bytes([50,100,10,10,1,1])
        self.assertEqual(list(gc_bins(data, 201, 200)), [(0,200,60,110), (200,201,1,1)])
        self.assertEqual(list(gc_bins(data, 201, 1000)), [(0,201,61,111)])
        # Off-grid estimates are explicitly proportional, including terminal tiles.
        self.assertEqual(list(gc_bins(bytes([100,100,0,100]), 200, 150)),
                         [(0,150,100.0,150.0), (150,200,0.0,50.0)])

    def test_truncated_and_trailing_count_payloads_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.bin.gz"
            for payload in (b"\x00", b"\x00\x01\x00"):
                with gzip_output(path) as stream:
                    write_count_header(stream, {"chr1":1})
                    stream.write(payload)
                with self.assertRaises(ValueError):
                    list(count_chromosomes(path))

    def test_cached_sources_must_match_the_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            path = cache / "fixture.dat"
            path.write_bytes(b"reference snapshot")
            metadata = {"url":"https://example.org/reference", "sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                        "bytes":path.stat().st_size, "retrieved_at_utc":"2026-10-24T00:00:00+00:00"}
            Path(str(path)+".http.json").write_text(json.dumps(metadata))
            args = Namespace(cache_dir=cache, download=False, refresh_sources=False)
            task = ("hg38", "fixture", path.name, metadata["url"], args, metadata)
            self.assertEqual(acquire(task)[2], metadata)
            path.write_bytes(b"unexpected update")
            with self.assertRaisesRegex(ValueError, "lock mismatch"):
                acquire(task)

    def test_gene_coverage_exact_with_nondivisor_and_terminal_bins(self):
        intervals = [(5, 20), (95, 105), (205, 207)]
        self.assertEqual(list(refseq_bins(intervals, 207, 100)),
                         [(0,100,0.2), (100,200,0.05), (200,207,2/7)])
        totals = list(refseq_bins(intervals, 207, 7))
        self.assertAlmostEqual(sum((e-s)*v for s,e,v in totals), 27)

    def test_refseq_deduplicates_isoforms_without_bridging_loci(self):
        def row(accession, start, end, name, chrom="chr1", strand="+"):
            return "\t".join(map(str, [0,accession,chrom,strand,start,end,start,end,1,
                                       f"{start},",f"{end},",0,name,"none","none","-1,"])) + "\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "refGene.gz"
            with gzip.open(path, "wt") as out:
                out.write(row("NM_1",10,30,"G1"))
                out.write(row("NM_2",20,40,"G1"))
                out.write(row("NR_3",50,60,"G1"))
                out.write(row("NM_4",35,55,"G2",strand="-"))
                out.write(row("NR_5",0,10,"G3",chrom="chr1_alt"))
            genes, coverage, stats = read_refseq(path, {"chr1":100})
            self.assertEqual(genes["chr1"], [(10,40,"G1","+"),(35,55,"G2","-"),(50,60,"G1","+")])
            self.assertEqual(coverage["chr1"], [(10,60)])
            self.assertEqual(stats["excluded_noncanonical_rows"], 1)
            self.assertEqual(stats["covered_bp"], 50)

    def test_export_cli_exactness_missing_values_and_hg38_default(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / "hg38"
            folder.mkdir()
            (folder / "chrom.sizes").write_text("chr1\t301\n")
            bases = "G"*50 + "A"*50 + "N"*100 + "C"*100 + "T"
            source = root / "source.2bit"
            fixture_2bit(source, bases)
            source_digest = hashlib.sha256(source.read_bytes()).hexdigest()
            (folder / "manifest.json").write_text(json.dumps({"sources":{"sequence":{"sha256":source_digest}}}))
            with gzip_output(folder / "gc.100bp.counts.bin.gz") as stream:
                write_count_header(stream, {"chr1":301})
                stream.write(bytes([50,100,0,0,100,100,0,1]))
            output = root / "gc.bedGraph"
            command = [sys.executable, str(SCRIPTS / "export_phasing_track.py"),
                       "--reference-dir", str(root), "--resolution", "100", "--output", str(output)]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_text().splitlines(), ["chr1\t0\t100\t0.5", "chr1\t100\t200\tnan",
                                                              "chr1\t200\t300\t1", "chr1\t300\t301\t0"])
            metadata = json.loads(Path(str(output)+".json").read_text())
            self.assertTrue(metadata["exact"])
            self.assertEqual(metadata["assembly"], "hg38")
            command[command.index("--resolution")+1] = "125"
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--two-bit", result.stderr)
            result = subprocess.run(command+["--allow-approximate"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(json.loads(Path(str(output)+".json").read_text())["exact"])
            result = subprocess.run(command+["--two-bit",str(source)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(json.loads(Path(str(output)+".json").read_text())["exact"])
            rows = [line.split() for line in output.read_text().splitlines()]
            self.assertTrue(math.isnan(float(rows[1][3])))  # 50 callable bases / 125 bp < 0.5
            self.assertAlmostEqual(float(rows[2][3]), 50/51)
            source.write_bytes(source.read_bytes()+b"different")
            result = subprocess.run(command+["--two-bit",str(source)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("checksum", result.stderr)


if __name__ == "__main__":
    unittest.main()
