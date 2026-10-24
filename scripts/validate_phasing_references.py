#!/usr/bin/env python3
"""Validate bundled hashes, all counts, interval coverage, and source agreement."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import random

from phasing_reference import (ASSEMBLIES, BIN_SIZE, ROOT, TwoBit,
                               count_chromosomes, merge_intervals, read_sizes)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def independently_decode(sequence, chrom, start, end):
    """Basewise decoder, independent of packed_gc and its bit-count algorithm."""
    length, nblocks, at = sequence.record(chrom)
    require(0 <= start < end <= length, "Invalid sequence query")
    bases = bytearray(b"TCAG"[(sequence.data[at+p//4] >> (6-2*(p%4))) & 3]
                      for p in range(start, end))
    for left, right in nblocks:
        a, b = max(start, left), min(end, right)
        if a < b:
            bases[a-start:b-start] = b"N" * (b-a)
    return bases.count(b"G") + bases.count(b"C"), len(bases) - bases.count(b"N")


# Translate each packed byte to the number of C/G bases it encodes. This supplies
# an independent genome-wide check of the builder's integer-mask/bit_count path.
GC_LOOKUP = bytes(sum(((b >> shift) & 3) in (1, 3) for shift in (6, 4, 2, 0))
                  for b in range(256))


def independent_raw_gc(sequence, at, start, end):
    first, last = (start + 3) // 4, end // 4
    if first > last:
        return sum(((sequence.data[at+p//4] >> (6-2*(p%4))) & 3) in (1, 3)
                   for p in range(start, end))
    total = 0
    for lo, hi in ((start, first*4), (last*4, end)):
        total += sum(((sequence.data[at+p//4] >> (6-2*(p%4))) & 3) in (1, 3)
                     for p in range(lo, hi))
    for offset in range(first, last, 1024 * 1024):
        chunk = sequence.data[at+offset:at+min(last, offset+1024*1024)].translate(GC_LOOKUP)
        total += sum(i * chunk.count(i) for i in range(1, 5))
    return total


def validate(assembly, root, source_cache=None):
    folder = root / assembly
    manifest = json.loads((folder / "manifest.json").read_text())
    require(manifest["assembly"] == assembly and manifest["gc_bin_size"] == BIN_SIZE, "Manifest mismatch")
    sizes = read_sizes(folder / "chrom.sizes")
    expected = [f"chr{i}" for i in range(1, 23 if assembly.startswith("hg") else 20)] + ["chrX", "chrY", "chrM"]
    require(list(sizes) == manifest["chromosomes"] == expected, "Unexpected chromosome set/order")
    lock = json.loads((root / "sources.lock.json").read_text())
    require(manifest["sources"] == lock[assembly], "Manifest/source lock disagreement")
    for name, info in manifest["files"].items():
        path = folder / name
        require(path.stat().st_size == info["bytes"] and digest(path) == info["sha256"], f"Checksum failure: {path}")
    if source_cache:
        for role, info in manifest["sources"].items():
            filename = {"sequence": f"{assembly}.2bit", "chrom_sizes": f"{assembly}.chrom.sizes",
                        "refseq": f"{assembly}.refGene.txt.gz"}[role]
            require(digest(source_cache / filename) == info["sha256"], f"Source checksum mismatch: {filename}")
    sequence = TwoBit(source_cache / f"{assembly}.2bit") if source_cache else None
    seen, total_bins = [], 0
    try:
        for chrom, length, data in count_chromosomes(folder / "gc.100bp.counts.bin.gz"):
            seen.append(chrom)
            require(sizes[chrom] == length, f"Length mismatch: {chrom}")
            bins = (length + BIN_SIZE - 1) // BIN_SIZE
            stats = manifest["gc_statistics"][chrom]
            gc_values, valid_values = data[::2], data[1::2]
            require(all(g <= v <= BIN_SIZE for g, v in zip(gc_values, valid_values)), f"Invalid counts: {chrom}")
            require(data[-2] <= data[-1] <= (length-1) % BIN_SIZE + 1, f"Terminal tile error: {chrom}")
            total_gc, total_valid = sum(gc_values), sum(valid_values)
            require(total_gc == stats["gc_bases"] and total_valid == stats["acgt_bases"], f"Count total mismatch: {chrom}")
            require(valid_values.count(0) == stats["no_sequence_bins"] and bins == stats["bins"], f"Bin total mismatch: {chrom}")
            if sequence:
                sequence_length, blocks, at = sequence.record(chrom)
                require(sequence_length == length, f"Source length mismatch: {chrom}")
                raw_gc = independent_raw_gc(sequence, at, 0, length)
                raw_gc -= sum(independent_raw_gc(sequence, at, a, b) for a, b in blocks)
                require(raw_gc == total_gc, f"Independent whole-chromosome GC mismatch: {chrom}")
                require(length - sum(b-a for a, b in blocks) == total_valid, f"Independent callable total mismatch: {chrom}")
                rng = random.Random(f"{assembly}/{chrom}")
                sample_tiles = {0, bins-1, *(rng.randrange(bins) for _ in range(30))}
                for a, b in blocks:
                    sample_tiles.update(t for t in (a//BIN_SIZE, b//BIN_SIZE) if t < bins)
                for tile in sorted(sample_tiles):
                    start, end = tile * BIN_SIZE, min(length, (tile+1)*BIN_SIZE)
                    expected_counts = independently_decode(sequence, chrom, start, end)
                    require(tuple(data[2*tile:2*tile+2]) == expected_counts, f"Independent tile mismatch: {chrom}:{start}")
            total_bins += bins
    finally:
        if sequence:
            sequence.close()
    require(seen == list(sizes), "GC chromosome order mismatch")
    loci = {chrom: [] for chrom in sizes}
    gene_count = 0
    rank = {chrom: i for i, chrom in enumerate(sizes)}
    previous = (-1, -1)
    with gzip.open(folder / "refseq.genes.bed.gz", "rt") as stream:
        for line in stream:
            chrom, start, end, symbol, score, strand = line.split()
            start, end = int(start), int(end)
            require(chrom in sizes and 0 <= start < end <= sizes[chrom], "Invalid gene locus")
            require(symbol and score == "0" and strand in ("+", "-"), "Invalid BED6 fields")
            key = rank[chrom], start
            require(key >= previous, "Unsorted genes")
            previous = key
            loci[chrom].append((start, end))
            gene_count += 1
    union = {chrom: merge_intervals(spans) for chrom, spans in loci.items()}
    cursors, active = {chrom: 0 for chrom in sizes}, {chrom: [] for chrom in sizes}
    previous_rank = -1
    with gzip.open(folder / "refseq.coverage.bedGraph.gz", "rt") as stream:
        for line in stream:
            chrom, start, end, value = line.split()
            start, end = int(start), int(end)
            require(chrom in sizes and rank[chrom] >= previous_rank, "Unsorted coverage chromosomes")
            previous_rank = rank[chrom]
            require(start == cursors[chrom] and start < end <= sizes[chrom], "Coverage gap, overlap, or bad boundary")
            require(value in ("0", "1"), "Nonbinary gene coverage")
            cursors[chrom] = end
            if value == "1":
                active[chrom].append((start, end))
    require(cursors == sizes, "Coverage does not reach all chromosome ends")
    require(active == union, "Coverage is not the union of BED gene loci")
    require(gene_count == manifest["refseq_statistics"]["gene_loci"], "Gene locus count mismatch")
    covered_bp = sum(b-a for spans in active.values() for a, b in spans)
    require(covered_bp == manifest["refseq_statistics"]["covered_bp"], "Gene coverage total mismatch")
    print(f"PASS {assembly}: {len(sizes)} chromosomes, {total_bins:,} exact 100 bp GC tiles, {gene_count:,} gene loci", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=ROOT)
    parser.add_argument("--assemblies", nargs="+", choices=ASSEMBLIES, default=list(ASSEMBLIES))
    parser.add_argument("--source-cache", type=Path)
    args = parser.parse_args()
    for assembly in args.assemblies:
        validate(assembly, args.reference_dir, args.source_cache)


if __name__ == "__main__":
    main()
