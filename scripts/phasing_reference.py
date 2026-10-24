#!/usr/bin/env python3
"""Standard-library readers and interval operations for bundled phasing data."""
from contextlib import contextmanager
import gzip
import json
import mmap
from pathlib import Path
import struct

MAGIC = b"GCCOUNT1"
BIN_SIZE = 100
ROOT = Path(__file__).resolve().parents[1] / "resources" / "compartment_phasing"
ASSEMBLIES = ("hg19", "hg38", "mm9", "mm10")


@contextmanager
def gzip_output(path, text=False):
    """Deterministic gzip: no timestamp or source filename in the header."""
    import io
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0, compresslevel=9) as zipped:
            if text:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="\n") as stream:
                    yield stream
            else:
                yield zipped


def read_sizes(path):
    result = {}
    with open(path) as stream:
        for line in stream:
            chrom, length = line.split()
            if chrom in result or int(length) <= 0:
                raise ValueError(f"Invalid chromosome sizes: {path}")
            result[chrom] = int(length)
    return result


def merge_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError("Invalid interval")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


class TwoBit:
    """Read UCSC version-0 2bit, either byte order; count GC excluding N blocks.

    Repeat-masked sequence is included. N blocks are explicitly removed from
    BOTH numerator and denominator, regardless of their underlying packed bits.
    """
    def __init__(self, path):
        self.file = open(path, "rb")
        self.data = mmap.mmap(self.file.fileno(), 0, access=mmap.ACCESS_READ)
        signature = self.data[:4]
        if signature == struct.pack("<I", 0x1A412743):
            self.order = "<"
        elif signature == struct.pack(">I", 0x1A412743):
            self.order = ">"
        else:
            self.close()
            raise ValueError("Not a UCSC 2bit file")
        version, count, reserved = struct.unpack_from(self.order + "III", self.data, 4)
        if version != 0 or reserved != 0:
            self.close()
            raise ValueError("Only UCSC version-0 2bit is supported")
        self.index = {}
        offset = 16
        for _ in range(count):
            size = self.data[offset]
            name = self.data[offset + 1:offset + 1 + size].decode("ascii")
            offset += 1 + size
            record, = struct.unpack_from(self.order + "I", self.data, offset)
            offset += 4
            if name in self.index:
                raise ValueError(f"Duplicate 2bit sequence: {name}")
            self.index[name] = record
        self.records = {}

    def close(self):
        self.data.close()
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def record(self, chrom):
        if chrom not in self.records:
            at = self.index[chrom]
            length, count = struct.unpack_from(self.order + "II", self.data, at)
            at += 8
            starts = struct.unpack_from(self.order + f"{count}I", self.data, at)
            at += count * 4
            sizes = struct.unpack_from(self.order + f"{count}I", self.data, at)
            at += count * 4
            mask_count, = struct.unpack_from(self.order + "I", self.data, at)
            at += 4 + mask_count * 8
            reserved, = struct.unpack_from(self.order + "I", self.data, at)
            at += 4
            if reserved or at + (length + 3) // 4 > len(self.data):
                raise ValueError(f"Invalid 2bit sequence record: {chrom}")
            blocks = [(s, s + n) for s, n in zip(starts, sizes)]
            if any(s < 0 or e > length or e <= s for s, e in blocks):
                raise ValueError(f"Invalid N blocks: {chrom}")
            if any(blocks[i][0] < blocks[i-1][1] for i in range(1, len(blocks))):
                raise ValueError(f"Unsorted/overlapping N blocks: {chrom}")
            self.records[chrom] = (length, blocks, at)
        return self.records[chrom]

    def packed_gc(self, packed_at, start, end):
        # A/C/G/T encoding has its low bit set exactly for C and G.
        left, right = start // 4, (end + 3) // 4
        count = 0
        if start % 4:
            stop = min(end, (left + 1) * 4)
            byte = self.data[packed_at + left]
            count += sum((byte >> (6 - 2 * (p % 4))) & 1 for p in range(start, stop))
            left += 1
        if end % 4 and left < right:
            byte = self.data[packed_at + right - 1]
            count += sum((byte >> (6 - 2 * p)) & 1 for p in range(end % 4))
            right -= 1
        if left < right:
            chunk = self.data[packed_at + left:packed_at + right]
            mask = (1 << (8 * len(chunk))) // 3  # 0x55 repeated
            count += (int.from_bytes(chunk, "big") & mask).bit_count()
        return count

    def counts(self, chrom, width=BIN_SIZE):
        length, blocks, packed_at = self.record(chrom)
        cursor = 0
        for start in range(0, length, width):
            end = min(length, start + width)
            gc, valid = self.packed_gc(packed_at, start, end), end - start
            while cursor < len(blocks) and blocks[cursor][1] <= start:
                cursor += 1
            j = cursor
            while j < len(blocks) and blocks[j][0] < end:
                a, b = max(start, blocks[j][0]), min(end, blocks[j][1])
                gc -= self.packed_gc(packed_at, a, b)
                valid -= b - a
                j += 1
            if not 0 <= gc <= valid <= end - start:
                raise ValueError(f"Invalid GC count: {chrom}:{start}-{end}")
            yield start, end, gc, valid


def write_count_header(stream, sizes):
    header = {"version": 1, "bin_size": BIN_SIZE, "columns": ["gc_bases", "acgt_bases"],
              "dtype": "uint8", "layout": "interleaved", "chromosomes": list(sizes.items())}
    encoded = json.dumps(header, separators=(",", ":")).encode("ascii")
    stream.write(MAGIC + struct.pack("<I", len(encoded)) + encoded)


def read_count_header(stream):
    if stream.read(8) != MAGIC:
        raise ValueError("Invalid GC-count magic")
    raw = stream.read(4)
    if len(raw) != 4:
        raise ValueError("Truncated GC-count header")
    length, = struct.unpack("<I", raw)
    if length > 1_000_000:
        raise ValueError("Invalid GC-count header size")
    header = json.loads(stream.read(length))
    if (header.get("version") != 1 or header.get("bin_size") != BIN_SIZE or
            header.get("columns") != ["gc_bases", "acgt_bases"] or
            header.get("dtype") != "uint8" or header.get("layout") != "interleaved"):
        raise ValueError("Unsupported GC-count format")
    return header


def count_chromosomes(path):
    """Yield one chromosome's interleaved counts at a time (about 5 MB)."""
    with gzip.open(path, "rb") as stream:
        header = read_count_header(stream)
        for chrom, length in header["chromosomes"]:
            count = (length + BIN_SIZE - 1) // BIN_SIZE
            data = stream.read(2 * count)
            if len(data) != 2 * count:
                raise ValueError(f"Truncated counts for {chrom}")
            yield chrom, length, data
        if stream.read(1):
            raise ValueError("Trailing GC-count data")
