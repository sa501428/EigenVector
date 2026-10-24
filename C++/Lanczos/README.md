**Lanczos Method with Selective Orthogonalization**  

This package provides tools for computing leading eigenvectors of contact matrices from Hi-C data.

## Building the Software

Use the [shared CMake build and regression tests](../README.md#build-and-test).
The Linux/macOS/MinGW wrapper scripts call the same CMake build. Install the
dependencies first; the scripts do not install system packages. `STRAW_DIR`
can override the path to the Straw C++ source directory.

## A/B orientation

All three programs use bundled GC content with a RefSeq agreement check for
hg19, hg38, mm9, and mm10. Suffixed genome IDs are recognized; chromosome names
and exact nuclear sizes confirm the ID or identify an inconclusive assembly.
`--assembly NAME` explicitly overrides detection, subject to size validation. GC is exact for 100 bp-aligned
resolutions; auto mode uses exact RefSeq coverage otherwise. `-p FILE` /
`--phasing-track FILE` supplies a custom activity bedGraph. All reference controls
are forwarded to both rescue phases. Unresolved/conflicting signs are left
unchanged and recorded in `.orientation.tsv` reports, including in quiet mode.
See [the shared phasing documentation](../README.md#ab-sign-orientation-at-any-resolution)
for options, confidence thresholds, output reports, and biological limitations.
Y and mitochondrial entries are retained, subject to adequate evidence.

## Usage

### Chromosome-parallel genome-wide analysis (`Lan.exe`)
```bash
./Lan.exe [options] <hicfile> <outbase> [resolution]

Options:
  --rescue     Select a high-resolution EV using a 50 kb reference
  -r <int>     Resolution in bp (default: 5000; can also be positional)
  -n <string>  Normalization method (default: SCALE)
  -o           Use observed instead of observed/expected (default: o/e)
  -j <int>     Concurrent chromosome processes (default: min(CPUs, 4))
  -T <int>     Threads per chromosome process (default: 1)
  -t <float>   Set tolerance (default: 1.0e-7)
  -e <float>   Set epsilon (default: 1.0e-8)
  -I <int>     Set maximum iterations (default: 200)
  -v <int>     Set verbosity level (default: 1)
  --coverage-z <float>
               High-resolution log1p(nonzero-count) z cutoff (default: -2.5)
  --no-coverage-filter
               Disable low-coverage row/column filtering
  --max-top1-energy <float>
               Reject localized rescue candidates above this fraction (default: 0.5)
```

Normal mode calculates EV1 independently for every chromosome and merges the
results into `<outbase>.wig`. Per-chromosome intermediates are removed after a
successful merge.
Missing parent directories in `<outbase>` are created automatically.
At unusually coarse resolutions, a chromosome may be too sparse to form an
eigensystem; that chromosome is zero-filled with a warning so the genome-wide
WIG remains structurally complete.

Before a requested-resolution calculation, the worker counts nonzero matrix
entries incident on every bin. It computes z-scores on `log1p(count)` among
nonzero bins, removes zero-count bins and bins below the default `-2.5` cutoff,
and removes every matrix entry touching a filtered bin. The 50 kb rescue
reference is intentionally left unfiltered. The log transform prevents the
strongly right-skewed count distribution from making the raw z-score cutoff
ineffective.

Example using all defaults (SCALE, o/e, 5 kb):

```bash
./Lan.exe sample.hic sample_lanczos
```

Example at 10 kb:

```bash
./Lan.exe sample.hic sample_lanczos 10000
```

#### Rescue mode

```bash
./Lan.exe --rescue sample.hic sample_lanczos_rescued 10000
```

For each chromosome, rescue mode calculates EV1 at 50 kb and EV1-EV10 at the
requested resolution. It expands each 50 kb call onto the requested-resolution
bins, then compares each high-resolution EV, in EV order, with that expanded
reference. The first vector whose
absolute Pearson correlation is greater than 0.8 is selected. The final A/B sign
requires adequate direct biological evidence; weak or conflicting evidence stays
unresolved. `applied_sign` and
`orientation_source` in the rescue TSV make this decision explicit. If no vector
passes, the best finite match is used and marked `fallback_best` in the report.
If a chromosome is too sparse to calculate a 50 kb reference (for example,
mitochondrial DNA), rescue retains high-resolution EV1 and reports
`no_low_resolution_reference`; its low-resolution reference section is zero-filled.

Before correlation-based selection, each high-resolution vector is checked for
localization. A vector is ineligible when its largest 1% of nonzero bins contain
more than 50% of its squared magnitude. If every finite candidate is localized,
the best match is retained only as an explicitly marked
`fallback_localized_no_eligible_vector` result.

Rescue mode produces:

- `<outbase>.wig`: rescued, genome-wide vector
- `<outbase>.lowres_50kb.wig`: genome-wide 50 kb reference
- `<outbase>.Ev1.wig` through `<outbase>.Ev10.wig`: unmodified high-resolution vectors
- `<outbase>.eigenvalues.tsv`: long-form per-chromosome eigenvalues and eigengap diagnostics
- `<outbase>.rescue.tsv`: selected EV, status, eigengap diagnostics, and all ten correlations

For the subspace containing the first `k` eigenvectors, the TSV diagnostics use

```text
gap_k       = lambda_k - lambda_(k+1)
gap_(k+1)   = lambda_(k+1) - lambda_(k+2)
min_ratio   = min(lambda_k / gap_k, lambda_(k+1) / gap_(k+1))
relative_error_estimate = tolerance * min_ratio
```

The two extra Ritz values maintained internally by the solver make these fields
available through `k=10`. The long-form eigenvalue TSV contains one row per
chromosome and `k`; the rescue TSV includes the same diagnostics for the
selected `k`. It also reports inverse participation ratio, top-1%-energy
fraction, 50 kb and 1 Mb spatial autocorrelation, localization flags, and the
coverage-filter threshold and retained/removed bin counts. The long-form
eigenvalue TSV includes these metrics for every `k`.

Use `--threshold` to change 0.8 and `--keep-temp` to retain per-chromosome
outputs and logs.

### Single-chromosome worker (`LanChr.exe`)

`LanChr.exe` is built for direct debugging and as the worker used by `Lan.exe`:

```bash
./LanChr.exe [options] <hicfile> <chromosome> <outbase> <resolution> [nv]
```

It defaults to SCALE, o/e, and one requested eigenvector. The solver retains up to two
additional internal Ritz vectors for diagnostics, but only `nv` vectors are
written. Exact rank exhaustion is handled without dividing by a zero residual;
unavailable extra eigenvalues are zero. When fewer independent
directions are available, the worker retries with fewer requested vectors.
Unavailable output vectors are zero-filled with a warning, and
`computed_eigenvectors` in the worker diagnostics records the available count.
Rescue reports `no_high_resolution_eigenvector` when no usable vector exists.
Eigenvalues, gaps, both conditioning ratios, and the tolerance-scaled
relative-error estimate are written to `<outbase>.eigenvalues.tsv`.

### Genome-Wide Analysis (GWev.exe)
```bash
./GWev.exe [options] <hicfile> <outbase> <resolution> [nv]

Options:
  -f           Use full matrix instead of inter-chromosomal only
  -t <float>   Set tolerance (default: 1.0e-7)
  -e <float>   Set epsilon (default: 1.0e-8)
  -I <int>     Set maximum iterations (default: 200)
  -T <int>     Set number of threads (default: 1)
  -v <int>     Set verbosity level (default: 1)
```

## Output Format
The programs generate eigenvector files in WIG format that can be visualized in genome browsers.
