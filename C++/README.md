# POSSUM C++ solvers

Power Method and Lanczos compute leading eigenvectors of the Pearson correlation
matrix of a symmetric contact matrix without materializing that correlation matrix.
Both methods now share sparse correlation statistics and A/B orientation logic.

## Build and test

Requires CMake 3.16+, a C++17 compiler, pthreads, curl, zlib, zstd, BLAS,
LAPACK and LAPACKE, plus the Straw C++ sources. Point `STRAW_DIR` at the directory
containing `straw.h` and `straw.cpp`:

```sh
cmake -S C++ -B C++/build -DCMAKE_BUILD_TYPE=Release -DSTRAW_DIR=/path/to/straw/C++
cmake --build C++/build --parallel
ctest --test-dir C++/build --output-on-failure
```

Run these commands from the repository root. CMake builds all six executables:
`ev.exe`, `GWevIntra.exe`, `mainEigen.exe`, `LanChr.exe`, `Lan.exe`, and `GWev.exe`.
It includes the additional Straw implementation files when present, including the
v10 reader and normalization-vector implementation. The platform scripts under
`Lanczos/` invoke this same build and place executables in `Lanczos/`.
They do not install system packages.

The test suite includes all four real phasing bundles, assembly aliases,
missing sequence, shortened terminal bins, off-grid resolutions, ambiguous evidence,
custom tracks, and global trans signs; comparisons
against independently constructed dense LAPACK Pearson eigenpairs (including
missing/constant bins and rank-one matrices); and complete production CLI runs
on a generated multi-chromosome, multi-resolution Hi-C file. Python integration
tests require Python 3.11+ and only the standard library.

## A/B sign orientation at any resolution

All five Hi-C entry points (`ev.exe`, `GWevIntra.exe`, `LanChr.exe`, `Lan.exe`,
`GWev.exe`) read the bundled GC/RefSeq references directly. Running a solver needs
no Python, sequence download, or pre-exported track. The sign convention is
**positive = A, negative = B**, conditional on adequate biological evidence.
The ASCII-only `mainEigen.exe` has no assembly/chromosome metadata and does not
perform biological sign orientation.

### Assembly and reference selection

The default is `--phasing-method auto`. Assembly detection reads the `.hic`
header's genome ID **and its complete chromosome names/sizes**, before chromosome
selection or solving:

1. Look for **hg19, hg38, mm9, or mm10 anywhere in the genome ID**, ignoring case.
   Prefixes and suffixes such as `project_HG38_no_alt` are supported. GRCh37/b37,
   GRCh38, NCBI37/NCBIM37, and GRCm38 aliases also support suffixes.
2. Confirm an unambiguous ID against every present canonical nuclear chromosome's
   exact size. If the ID is missing, unrecognized, ambiguous, or inconsistent,
   compare that chromosome metadata against all four bundles and use a **unique**
   matching assembly. A contradictory ID corrected by sizes emits a diagnostic.
3. Record the selected assembly and detection source. If no unique size match
   exists for an unknown ID, leave orientation unresolved. A recognizable ID with
   contradictory sizes and no unique fallback is an error. An explicit
   `--assembly NAME` override must also match the nuclear chromosome sizes.

Chromosome aliases with/without `chr` are normalized. Alternate, random, and
unplaced scaffolds are ignored for identification. Mitochondrial sizes do not
identify the nuclear assembly, since mitochondrial references can differ within
one nuclear build (notably b37 versus UCSC hg19). A mismatched mitochondrial
reference gets `chromosome_length_mismatch` and remains unresolved, while the
nuclear chromosomes can still be phased. A custom activity track bypasses bundled
assembly/size requirements because its coordinate system is user-supplied.

The export helper still defaults to **hg38**. For direct library callers without
any chromosome metadata and with an empty ID, hg38 is also the default. Hi-C CLI
runs always inspect the stored chromosome metadata and can therefore recover
hg19/mm9/mm10 even when the ID is missing or has a custom label.

References are found through `--phasing-dir DIR`, then `POSSUM_PHASING_DIR`, then
the path configured by CMake (the checkout's `resources/compartment_phasing` by
default). The directory must contain the assembly subdirectories. Relocated
binaries need the reference directory too. Manual builds without the CMake path
can set the environment variable or pass the directory explicitly. Missing or
malformed required data are errors; they do not silently disable phasing.

The [bundled reference documentation](../resources/compartment_phasing/README.md)
describes source snapshots, checksums, storage format, rebuilding, and validation.
The historical 100 kb A anchors are no longer used by the orientation helper.

### Default decision

1. Sum GC and valid A/C/G/T base counts into the actual Hi-C bins. Unknown sequence
   is excluded from both counts. GC is missing when fewer than 50% of bases are
   valid A/C/G/T; repeat-masked sequence remains valid.
2. Calculate the exact fraction of each bin covered by the union of RefSeq gene
   bodies. Overlapping transcripts and genes do not count their bases twice.
3. Calculate Spearman correlations using paired finite reference/EV bins. Zero EV
   values are excluded because these solvers use zero for masked/constant contact
   rows. Pearson correlations are also recorded for inspection.
4. Use GC when it passes the evidence thresholds. If GC is weak or unavailable,
   use RefSeq when it passes. If both pass but disagree in sign, mark the result
   `conflicting_references` and leave the vector unchanged.

By default a reference needs **20 paired bins** and **absolute Spearman correlation
at least 0.1**. Constant tracks, zero correlations, and insufficient evidence remain
unresolved. These are configurable engineering thresholds, not calibrated
statistical significance tests: spatially adjacent bins are correlated. Set
`--phasing-min-bins N` (at least 2) and `--phasing-min-correlation X` (0 through 1)
when justified for a dataset. Reducing them can admit unstable orientations.

GC is exact when the resolution is divisible by **100 bp**. At other resolutions,
auto mode uses exact RefSeq coverage without approximating GC. Forced
`--phasing-method gc` rejects such requests unless
`--phasing-allow-approximate-gc` is explicitly given; that option estimates the
partial boundary tiles by overlap and records `gc_exact=0`. An exact custom GC
track calculated from sequence is another option. RefSeq coverage is exact for
any positive integer resolution supported by the Hi-C map.

Use `--phasing-method gc` or `refseq` to select one reference alone, or `none` to
preserve solver signs without needing bundled data. Fine reference tiles do not
establish reliable compartment inference at fine resolution: contact coverage and
the identity of the component still need inspection.

### Custom activity tracks

`-p FILE` / `--phasing-track FILE` overrides the bundled references in auto mode,
including for unsupported assemblies and noncanonical chromosomes. It is rejected
in combination with an explicit `gc`, `refseq`, or `none` method. Use a plain
four-column bedGraph with zero-based, half-open coordinates on the same assembly:

```text
# chromosome start end activity
chr1 0 1000 2.4
chr1 1000 2000 0.3
chrY 0 1000 1.8
```

Higher values must correspond to more active chromatin. Intervals can have any
width, are sorted automatically, and may not overlap. Values are averaged using
interval overlap lengths, then oriented by Spearman correlation under the same
evidence thresholds. Missing/nonfinite values are excluded, finite zero activity
is valid, and absent chromosomes do not fall back to bundled data. A sparse list
of gene intervals with only value 1 needs explicit intergenic zeros to represent
gene density. Chromosome aliases `1`/`chr1`, `X`/`chrX`, `Y`/`chrY`, and
`M`/`MT`/`chrM`/`chrMT` are matched consistently.

### Diagnostics and genome-wide behavior

Every Hi-C output gets an **orientation TSV** with assembly and its origin
(`map_header`, `chromosome_sizes`, or `override`; with unresolved/default origins
when applicable), chromosome, resolution, component,
status, reference source/path, pre-flip Spearman/Pearson correlations, paired-bin
counts, evidence thresholds, and GC exactness. `applied_sign` is the actual factor
(+1 or -1); `decision_sign=0` identifies an unresolved or disabled decision, whose
vector is unchanged. An unresolved status must not be interpreted as a confirmed
A/B assignment. Eigenvector text uses 17 significant digits so saved values retain
the precision used in correlation diagnostics.

| Program | Orientation report |
| --- | --- |
| `ev.exe`, `GWevIntra.exe` | `<outfile>.orientation.tsv` |
| `LanChr.exe`, `GWev.exe`, `Lan.exe` | `<outbase>.orientation.tsv` |
| `Lan.exe --rescue` | Also `<outbase>.candidates.orientation.tsv` and `<outbase>.lowres_50kb.orientation.tsv` |

Cis vectors are phased per chromosome. Trans (`GWev.exe`) pools reference/EV pairs
across chromosomes to calculate one correlation and applies one **global** sign
per vector, preserving relative chromosome signs. Y and mitochondrial entries are
retained; this does not establish that conventional A/B interpretation is valid
for every chromosome or cell state.

Rescue retains the existing 50 kb-guided component selection and localization
checks. It uses biological reference evidence for the selected vector's final
sign; a 50 kb match alone does not resolve weak or conflicting A/B evidence. The
final orientation TSV is relative to the selected worker vector; the candidate
TSV preserves the workers' original phasing decisions. Both phases receive the
same assembly, reference, and threshold controls, and reports survive cleanup.

```sh
# Map metadata selects the bundled assembly automatically.
C++/build/ev.exe sample.hic chr1 chr1.ev.txt 500
C++/build/GWevIntra.exe sample.hic power.wig 1000
C++/build/LanChr.exe sample.hic chr1 chr1 500
C++/build/Lan.exe sample.hic lanczos 1000

# Explicit assembly or reference control.
C++/build/Lan.exe --assembly hg38 --rescue sample.hic rescued 500
C++/build/Lan.exe --phasing-method refseq sample.hic gene_phased 7777
C++/build/Lan.exe --phasing-track activity.bedgraph sample.hic custom_phased 500
```

The requested resolution must be present in the Hi-C file. See the CLI `-h`
output for normalization and solver settings, which are independent of phasing.
