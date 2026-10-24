#!/usr/bin/env python3
"""Build bundled references from checksummed UCSC sequence and RefSeq tables.

    python3 scripts/build_phasing_references.py --download --cache-dir /tmp/phasing

Large downloaded source files stay in the external cache, outside the repository.
Rebuilds enforce sources.lock.json; --refresh-sources explicitly updates the lock.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

from phasing_reference import (ASSEMBLIES, BIN_SIZE, ROOT, TwoBit, gzip_output,
                               merge_intervals, read_sizes, write_count_header)

HOST = "https://hgdownload.soe.ucsc.edu/goldenPath"


def sha256(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def source_specs(assembly):
    return {"sequence": (f"{assembly}.2bit", f"{HOST}/{assembly}/bigZips/{assembly}.2bit"),
            "chrom_sizes": (f"{assembly}.chrom.sizes", f"{HOST}/{assembly}/bigZips/{assembly}.chrom.sizes"),
            "refseq": (f"{assembly}.refGene.txt.gz", f"{HOST}/{assembly}/database/refGene.txt.gz")}


def acquire(task):
    assembly, role, name, url, args, locked = task
    path = args.cache_dir / name
    metadata_path = path.with_name(path.name + ".http.json")
    if not path.exists():
        if not args.download:
            raise ValueError(f"Missing {path}; run with --download")
        print(f"Downloading {url}", flush=True)
        temporary = path.with_name(path.name + ".part")
        request = urllib.request.Request(url, headers={"User-Agent": "possum-phasing-reference-builder/1"})
        with urllib.request.urlopen(request, timeout=120) as response, open(temporary, "wb") as out:
            expected_length = response.headers.get("Content-Length")
            size = 0
            while chunk := response.read(4 * 1024 * 1024):
                out.write(chunk)
                size += len(chunk)
            if expected_length is not None and size != int(expected_length):
                raise ValueError(f"Incomplete download: {url}")
            metadata = {"url": url, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                        "last_modified": response.headers.get("Last-Modified"),
                        "etag": response.headers.get("ETag")}
        temporary.replace(path)
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    if not metadata_path.exists():
        raise ValueError(f"Source provenance missing: {metadata_path}; re-download into a fresh cache")
    metadata = json.loads(metadata_path.read_text())
    if metadata["url"] != url:
        raise ValueError(f"Source URL mismatch: {path}")
    metadata.update({"sha256": sha256(path), "bytes": path.stat().st_size})
    if locked and not args.refresh_sources:
        for key in ("url", "sha256", "bytes"):
            if metadata[key] != locked[key]:
                raise ValueError(f"Source lock mismatch for {assembly}/{role}: {key}")
        # Preserve recorded retrieval information on byte-identical rebuilds.
        metadata = locked
    print(f"Verified {assembly}/{role}: {metadata['bytes']:,} bytes", flush=True)
    return assembly, role, metadata


def read_refseq(path, sizes):
    """Union transcript spans per symbol/chromosome/strand/connected locus.

    This never bridges disjoint copies of a symbol across a genomic gap. Coverage
    is then a union across ALL loci, so overlapping genes/isoforms count once.
    """
    grouped, accepted, skipped, coding, noncoding = {}, 0, 0, 0, 0
    with gzip.open(path, "rt") as stream:
        for number, line in enumerate(stream, 1):
            row = line.rstrip("\n").split("\t")
            if len(row) != 16:
                raise ValueError(f"Invalid refGene row {number}: expected 16 columns")
            accession, chrom, strand = row[1:4]
            if chrom not in sizes:
                skipped += 1
                continue
            start, end = int(row[4]), int(row[5])
            if strand not in ("+", "-") or not 0 <= start < end <= sizes[chrom]:
                raise ValueError(f"Invalid refGene interval at row {number}")
            if not accession.startswith(("NM_", "NR_")):
                raise ValueError(f"Unexpected refGene accession at row {number}: {accession}")
            name = row[12] or accession
            grouped.setdefault((chrom, name, strand), []).append((start, end))
            accepted += 1
            coding += int(int(row[6]) != int(row[7]))
            noncoding += int(int(row[6]) == int(row[7]))
    genes, covered = {c: [] for c in sizes}, {c: [] for c in sizes}
    for (chrom, name, strand), spans in grouped.items():
        for start, end in merge_intervals(spans):
            genes[chrom].append((start, end, name, strand))
            covered[chrom].append((start, end))
    for chrom in sizes:
        genes[chrom].sort()
        covered[chrom] = merge_intervals(covered[chrom])
        if chrom != "chrM" and not covered[chrom]:
            raise ValueError(f"No RefSeq genes on expected chromosome {chrom}")
    return genes, covered, {"accepted_transcript_rows": accepted, "excluded_noncanonical_rows": skipped,
                            "coding_transcript_rows": coding, "noncoding_transcript_rows": noncoding,
                            "gene_loci": sum(map(len, genes.values())),
                            "covered_bp": sum(e-s for spans in covered.values() for s, e in spans)}


def build(assembly, args, sources):
    out = args.output_dir / assembly
    out.mkdir(parents=True, exist_ok=True)
    all_sizes = read_sizes(args.cache_dir / f"{assembly}.chrom.sizes")
    names = [f"chr{i}" for i in range(1, 23 if assembly.startswith("hg") else 20)] + ["chrX", "chrY", "chrM"]
    sizes = {c: all_sizes[c] for c in names}
    (out / "chrom.sizes").write_text("".join(f"{c}\t{n}\n" for c, n in sizes.items()))
    genes, covered, refseq_stats = read_refseq(args.cache_dir / f"{assembly}.refGene.txt.gz", sizes)
    with gzip_output(out / "refseq.genes.bed.gz", text=True) as stream:
        for chrom in sizes:
            for start, end, symbol, strand in genes[chrom]:
                stream.write(f"{chrom}\t{start}\t{end}\t{symbol}\t0\t{strand}\n")
    with gzip_output(out / "refseq.coverage.bedGraph.gz", text=True) as stream:
        for chrom, length in sizes.items():
            cursor = 0
            for start, end in covered[chrom]:
                if start > cursor:
                    stream.write(f"{chrom}\t{cursor}\t{start}\t0\n")
                stream.write(f"{chrom}\t{start}\t{end}\t1\n")
                cursor = end
            if cursor < length:
                stream.write(f"{chrom}\t{cursor}\t{length}\t0\n")
    gc_stats = {}
    with TwoBit(args.cache_dir / f"{assembly}.2bit") as sequence:
        with gzip_output(out / "gc.100bp.counts.bin.gz") as stream:
            write_count_header(stream, sizes)
            for chrom, length in sizes.items():
                if sequence.record(chrom)[0] != length:
                    raise ValueError(f"2bit/chrom.sizes mismatch: {assembly}/{chrom}")
                buffer, total_gc, total_valid, empty, bins = bytearray(), 0, 0, 0, 0
                for start, end, gc, valid in sequence.counts(chrom):
                    buffer.extend((gc, valid))
                    total_gc += gc
                    total_valid += valid
                    bins += 1
                    empty += valid == 0
                    if len(buffer) >= 1024 * 1024:
                        stream.write(buffer)
                        buffer.clear()
                stream.write(buffer)
                gc_stats[chrom] = {"length": length, "bins": bins, "gc_bases": total_gc,
                                   "acgt_bases": total_valid, "no_sequence_bins": empty}
                print(f"Built {assembly}/{chrom}: {bins:,} GC bins", flush=True)
    files = {}
    for path in sorted(out.iterdir()):
        if path.name == "manifest.json":
            continue
        files[path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    manifest = {"schema_version": 1, "assembly": assembly, "gc_bin_size": BIN_SIZE,
                "coordinates": "zero-based, half-open", "chromosomes": list(sizes),
                "sequence_policy": "canonical chromosomes, including X/Y/M; soft masks included; N excluded",
                "refseq_policy": "UCSC refGene NM_/NR_ transcripts; connected gene loci; union gene-body coverage",
                "sources": sources, "gc_statistics": gc_stats, "refseq_statistics": refseq_stats,
                "files": files}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Finished {assembly}: {sum(v['bytes'] for v in files.values()):,} bundled bytes", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", nargs="+", choices=ASSEMBLIES, default=list(ASSEMBLIES))
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--refresh-sources", action="store_true", help="Explicitly accept new source checksums")
    args = parser.parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    lock_path = args.output_dir / "sources.lock.json"
    lock = json.loads(lock_path.read_text()) if lock_path.exists() else {}
    tasks = [(a, role, name, url, args, lock.get(a, {}).get(role))
             for a in args.assemblies for role, (name, url) in source_specs(a).items()]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for assembly, role, metadata in pool.map(acquire, tasks):
            lock.setdefault(assembly, {})[role] = metadata
    lock_path.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    for assembly in args.assemblies:
        build(assembly, args, lock[assembly])


if __name__ == "__main__":
    main()
