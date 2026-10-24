#!/usr/bin/env python3
"""Export a bundled assembly reference as a standard four-column bedGraph.

GC aggregation is exact for resolutions divisible by 100 bp. Other resolutions
require --two-bit for exact counts, or explicit --allow-approximate for uniform
overlap estimates within boundary tiles. RefSeq coverage is exact at any size.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

from phasing_reference import (ASSEMBLIES, BIN_SIZE, ROOT, TwoBit,
                               count_chromosomes, read_sizes)


def gc_bins(data, length, resolution):
    """Overlap-weighted integer counts; integral whenever boundaries align."""
    for start in range(0, length, resolution):
        end = min(start + resolution, length)
        gc, valid = 0, 0
        first = start // BIN_SIZE
        last = (end - 1) // BIN_SIZE
        if start % BIN_SIZE == 0 and (end % BIN_SIZE == 0 or end == length):
            # Fast exact path: integer sums, not mean tile percentages.
            yield start, end, sum(data[2*first:2*(last+1):2]), sum(data[2*first+1:2*(last+1):2])
            continue
        for tile in range(first, last + 1):
            left, right = tile * BIN_SIZE, min((tile + 1) * BIN_SIZE, length)
            overlap = min(end, right) - max(start, left)
            if overlap == right - left:
                gc += data[2 * tile]
                valid += data[2 * tile + 1]
            else:
                weight = overlap / (right - left)
                gc += data[2 * tile] * weight
                valid += data[2 * tile + 1] * weight
        yield start, end, gc, valid


def refseq_bins(intervals, length, resolution):
    cursor = 0
    for start in range(0, length, resolution):
        end = min(start + resolution, length)
        while cursor < len(intervals) and intervals[cursor][1] <= start:
            cursor += 1
        j, covered = cursor, 0
        while j < len(intervals) and intervals[j][0] < end:
            a, b = intervals[j]
            covered += min(b, end) - max(a, start)
            j += 1
        yield start, end, covered / (end - start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", choices=ASSEMBLIES, default="hg38")
    parser.add_argument("--reference-dir", type=Path, default=ROOT)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--track", choices=("gc", "refseq", "acgt-fraction"), default="gc")
    parser.add_argument("--chromosomes", nargs="+")
    parser.add_argument("--output", type=Path, required=True, help=".gz suffix enables gzip")
    parser.add_argument("--min-acgt-fraction", type=float, default=0.5,
                        help="GC is NaN below this sequence coverage (default: 0.5)")
    parser.add_argument("--allow-approximate", action="store_true")
    parser.add_argument("--two-bit", type=Path, help="Original locked assembly .2bit for exact arbitrary bins")
    args = parser.parse_args()
    if args.resolution <= 0 or not 0 <= args.min_acgt_fraction <= 1:
        parser.error("resolution must be positive and min-acgt-fraction must be in [0,1]")
    folder = args.reference_dir / args.assembly
    sizes = read_sizes(folder / "chrom.sizes")
    selected = set(args.chromosomes or sizes)
    if not selected or selected - sizes.keys():
        parser.error(f"Unknown chromosomes: {sorted(selected - sizes.keys())}")
    sequence_track = args.track in ("gc", "acgt-fraction")
    approximate = sequence_track and args.resolution % BIN_SIZE != 0 and not args.two_bit
    if approximate and not args.allow_approximate:
        parser.error("GC bins must align to 100 bp; use --two-bit for exact values or --allow-approximate")
    manifest = json.loads((folder / "manifest.json").read_text())
    if args.two_bit and sequence_track:
        with open(args.two_bit, "rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != manifest["sources"]["sequence"]["sha256"]:
            parser.error("--two-bit does not match the locked assembly sequence checksum")
    print(f"Exporting {args.assembly} {args.track}: {'approximate boundary tiles' if approximate else 'exact'}", file=sys.stderr)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if args.output.suffix == ".gz" else open
    rows = 0
    with opener(args.output, "wt") as out:
        def emit(chrom, start, end, value):
            nonlocal rows
            out.write(f"{chrom}\t{start}\t{end}\t{value:.10g}\n")
            rows += 1

        def emit_gc(chrom, records):
            for start, end, gc, valid in records:
                if args.track == "acgt-fraction":
                    value = valid / (end - start)
                else:
                    value = gc / valid if valid > 0 and valid / (end-start) >= args.min_acgt_fraction else float("nan")
                emit(chrom, start, end, value)

        if sequence_track and args.two_bit:
            with TwoBit(args.two_bit) as sequence:
                for chrom, length in sizes.items():
                    if chrom in selected:
                        if sequence.record(chrom)[0] != length:
                            raise ValueError(f"Sequence length mismatch: {chrom}")
                        emit_gc(chrom, sequence.counts(chrom, args.resolution))
        elif sequence_track:
            for chrom, length, data in count_chromosomes(folder / "gc.100bp.counts.bin.gz"):
                if length != sizes[chrom]:
                    raise ValueError(f"GC length mismatch: {chrom}")
                if chrom in selected:
                    emit_gc(chrom, gc_bins(data, length, args.resolution))
        else:
            intervals = {chrom: [] for chrom in sizes}
            with gzip.open(folder / "refseq.coverage.bedGraph.gz", "rt") as track:
                for line in track:
                    chrom, start, end, value = line.split()
                    if value == "1":
                        intervals[chrom].append((int(start), int(end)))
            for chrom, length in sizes.items():
                if chrom in selected:
                    for start, end, value in refseq_bins(intervals[chrom], length, args.resolution):
                        emit(chrom, start, end, value)
    metadata = {"assembly": args.assembly, "track": args.track, "resolution": args.resolution,
                "exact": not approximate, "chromosomes": [c for c in sizes if c in selected],
                "min_acgt_fraction": args.min_acgt_fraction if args.track == "gc" else None,
                "rows": rows, "source_manifest_sha256": hashlib.sha256((folder / "manifest.json").read_bytes()).hexdigest()}
    with open(args.output, "rb") as stream:
        metadata["output_sha256"] = hashlib.file_digest(stream, "sha256").hexdigest()
    args.output.with_name(args.output.name + ".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
