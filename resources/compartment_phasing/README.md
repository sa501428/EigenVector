# Bundled compartment phasing references

These reference data are for orienting Hi-C compartment eigenvectors with higher
activity-associated reference values on the positive (A) side. They do not contain
preassigned A/B compartments. All Hi-C solvers now read these bundles directly,
using GC with a RefSeq agreement check by default. See the
[solver options and diagnostics](../../C++/README.md#ab-sign-orientation-at-any-resolution).

The four supported assemblies are **hg19, hg38, mm9, and mm10**. The export helper
defaults to **hg38**. The solvers recognize these assembly names inside suffixed
genome IDs and confirm them against the stored chromosome sizes. Missing,
ambiguous, unrecognized, or inconsistent IDs fall back to a unique match of the
canonical nuclear chromosome names/sizes against these four bundles. An explicit
`--assembly` override is also size-checked. Unresolved signatures do not receive
an assumed biological A/B orientation. See the solver documentation for
mitochondrial-reference handling and custom tracks.

## Contents per assembly

| File | Content |
| --- | --- |
| `chrom.sizes` | Canonical chromosomes in numeric order, then X, Y, M. |
| `gc.100bp.counts.bin.gz` | Exact counts of G+C and A+C+G+T in every 100 bp tile. |
| `refseq.genes.bed.gz` | BED6 gene loci derived from RefSeq transcript spans. |
| `refseq.coverage.bedGraph.gz` | Complete 0/1 gene-body union coverage, including explicit zero intervals. |
| `manifest.json` | Source URLs, retrieval dates, source/output SHA-256 checksums, policies, and summary counts. |

`sources.lock.json` pins every upstream input. No full genome sequence or raw
annotation download is committed. Gzip is standard gzip, not BGZF/tabix-indexed.
All genomic coordinates are **zero-based, half-open**. Noncanonical scaffolds,
alternate loci, patches, and random contigs are excluded. X/Y/M are present for
completeness; inclusion does not assert that conventional A/B analysis is suitable
for every chromosome or experimental context.

The GC count file is deliberately small and self-describing. After gzip
decompression its layout is: 8 ASCII bytes `GCCOUNT1`, a 4-byte little-endian JSON
header length, that many UTF-8 JSON header bytes, then interleaved unsigned-byte
`gc_bases, acgt_bases` pairs. The header lists the chromosome names and lengths in
storage order. Each chromosome has `ceil(length / 100)` pairs with no padding;
the last pair describes the shortened terminal tile. Counts never exceed 100, so
one byte per count is lossless. There are no floating-point values to round or
average. `scripts/phasing_reference.py` provides a streaming reader, and the
export helper converts this storage format to standard bedGraph for other tools.

## Use at a requested resolution

Only Python 3.11+ and its standard library are required. Run from the repository
root (or use an absolute script path):

```sh
# The assembly defaults to hg38. Choose the actual Hi-C resolution.
python3 scripts/export_phasing_track.py --resolution 5000 --track gc \
  --output /tmp/hg38.gc.5kb.bedGraph

python3 scripts/export_phasing_track.py --assembly mm10 --resolution 500 \
  --track refseq --output /tmp/mm10.refseq.500bp.bedGraph

# Optional companion track to inspect reference sequence coverage.
python3 scripts/export_phasing_track.py --assembly hg19 --resolution 10000 \
  --track acgt-fraction --output /tmp/hg19.acgt.10kb.bedGraph
```

The output is standard four-column `chrom start end value` bedGraph. Use a `.gz`
suffix for compressed output. A JSON sidecar records the assembly, resolution,
exactness, reference manifest checksum, output checksum, and filtering settings.
Restrict output with, for example, `--chromosomes chr1 chr2`.

GC is calculated as **sum(G+C) / sum(A+C+G+T)**. It is not the unweighted mean of
tile percentages. Repeat-masked A/C/G/T bases are retained. Unknown N bases are
excluded from both numerator and denominator. All-N bins are `nan`, never zero GC.
By default, GC bins with less than 50% A/C/G/T sequence are also `nan`; change this
with `--min-acgt-fraction`. Gene coverage is the number of bases in the union of
gene bodies divided by the actual bin width, including a shortened terminal bin.
It measures annotation coverage, not expression or the number of active genes.

## Resolution guarantees

**GC is exact at any resolution that is a multiple of 100 bp**, with bins starting
at zero. This includes 100 bp, 500 bp, 1 kb, 2.5 kb, 5 kb, 7 kb, 10 kb, 25 kb,
and standard larger Hi-C resolutions. Chromosome terminal tiles retain their true
length. The counts are lossless at 100 bp, but do not preserve base order inside a
tile. They cannot provide exact GC for arbitrary boundaries or bins below 100 bp.

For exact GC at other resolutions, supply the original checksummed `.2bit` file:

```sh
python3 scripts/export_phasing_track.py --assembly hg38 --resolution 125 \
  --two-bit /tmp/possum-phasing-sources/hg38.2bit \
  --output /tmp/hg38.gc.125bp.bedGraph
```

Alternatively, `--allow-approximate` explicitly enables overlap-proportional
estimates within the boundary tiles; the sidecar marks them as approximate. The
helper refuses such GC requests by default. RefSeq coverage uses exact intervals
and is exact at **any positive integer resolution**, without sequence downloads.
The data's fine resolution does not by itself establish reliable compartment
inference at that resolution: Hi-C coverage and biological signal still matter.

## RefSeq annotation policy

The source for all four assemblies is UCSC's `refGene` database table: RefSeq
transcripts aligned by UCSC, **not** the distinct NCBI-aligned `ncbiRefSeq` track.
This gives a consistent source for mm9 as well as the other three builds. Included
transcripts have `NM_` or `NR_` accessions (coding and noncoding); the source tables
and SHA-256 hashes identify the downloaded snapshot, not a claimed numbered NCBI
annotation release. All isoforms are included in coverage without counting their
overlapping bases more than once.

For BED6 loci, spans are grouped by gene symbol, chromosome, and strand, then only
overlapping/touching spans are merged. Disjoint loci with the same symbol remain
separate. The coverage track is a further union across all gene loci, irrespective
of strand. Introns are included; intergenic regions have explicit zero coverage.
These loci are for coverage and display; BED row counts should not be interpreted
as a definitive count of distinct biological genes.

## Rebuild and verification

```sh
# Downloads about 3 GB of sequence plus small annotations into an external cache.
python3 scripts/build_phasing_references.py --download \
  --cache-dir /tmp/possum-phasing-sources

# Offline rebuild with cached inputs; enforces the recorded source checksums.
python3 scripts/build_phasing_references.py \
  --cache-dir /tmp/possum-phasing-sources

python3 scripts/validate_phasing_references.py
python3 -m unittest discover -s tests -p 'test_phasing_references.py'
```

The builder has `--assemblies` and `--output-dir` options. An upstream checksum
change causes failure; use a fresh cache with `--download --refresh-sources` only
when intentionally replacing the reference snapshot. When used, review the changed
source lock and manifests. Existing cached inputs are reused rather than fetched
again, including when `--refresh-sources` is passed.
Gzip headers are deterministic, so identical inputs on the same Python/zlib
runtime produce identical output bytes. The raw .2bit downloads stay outside Git
and may be deleted after building; retain them if exact off-grid GC export is
needed. Rebuilds use no shell-downloaded executables or third-party Python code.

The validator checks output hashes, every GC count and genomic boundary, chromosome
totals, complete nonoverlapping coverage, and equivalence of the BED loci with the
coverage union. Supply `--source-cache /tmp/possum-phasing-sources` to additionally
compare counts to independent sequence decoding and verify source checksums.

## Sources and method context

- [UCSC 2bit and BED format specifications](https://genome.ucsc.edu/FAQ/FAQformat.html).
- [UCSC gene annotation FAQ](https://genome.ucsc.edu/FAQ/FAQgenes.html), including
  the distinction between UCSC RefSeq and NCBI RefSeq and the hg19 caveat.
- [UCSC downloads](https://hgdownload.soe.ucsc.edu/downloads.html).
- [Cooltools compartments tutorial](https://cooltools.readthedocs.io/en/latest/notebooks/compartments_and_saddles.html)
  for GC phasing; [HiCExplorer hicPCA](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicPCA.html)
  for gene/experimental-track phasing.

These are derived biological reference data. UCSC describes GenBank/RefSeq data
distribution in its assembly download README files; source data are not represented
as original software authored by this project. See the upstream
[conditions of use](https://genome.ucsc.edu/conditions.html).
