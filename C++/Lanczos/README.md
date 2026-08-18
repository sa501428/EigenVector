**Lanczos Method with Selective Orthogonalization**  

This package provides tools for computing leading eigenvectors of contact matrices from Hi-C data.

## Building the Software

### Prerequisites

The software requires:
- C++ compiler (GCC/G++ 4.8 or later)
- OpenBLAS
- LAPACK/LAPACKE
- libcurl
- zlib
- straw library (place in `~/straw` or modify build script)

### Platform-Specific Build Instructions

#### Linux
```bash
# Make the script executable
chmod +x build_linux.sh
# Run the build script
./build_linux.sh
```

The script will automatically install required dependencies using apt-get.

#### macOS
```bash
# Make the script executable
chmod +x build_mac.sh
# Run the build script
./build_mac.sh
```

The script will use Homebrew to install required dependencies.

#### Windows
1. Install MinGW-w64 from [WinLibs](https://winlibs.com/)
2. Add MinGW-w64 bin directory to your PATH
3. Run the build script:
```cmd
build_windows.bat
```

You may need to modify the paths in `build_windows.bat` to match your MinGW-w64 installation.

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
```

Normal mode calculates EV1 independently for every chromosome and merges the
results into `<outbase>.wig`. Per-chromosome intermediates are removed after a
successful merge.
Missing parent directories in `<outbase>` are created automatically.
At unusually coarse resolutions, a chromosome may be too sparse to form an
eigensystem; that chromosome is zero-filled with a warning so the genome-wide
WIG remains structurally complete.

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
absolute Pearson correlation is greater than 0.8 is selected; a negative match is
sign-aligned to the 50 kb reference. If no vector passes, the best finite match
is used and marked `fallback_best` in the report.
If a chromosome is too sparse to calculate a 50 kb reference (for example,
mitochondrial DNA), rescue retains high-resolution EV1 and reports
`no_low_resolution_reference`; its low-resolution reference section is zero-filled.

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
selected `k`.

Use `--threshold` to change 0.8 and `--keep-temp` to retain per-chromosome
outputs and logs.

### Single-chromosome worker (`LanChr.exe`)

`LanChr.exe` is built for direct debugging and as the worker used by `Lan.exe`:

```bash
./LanChr.exe [options] <hicfile> <chromosome> <outbase> <resolution> [nv]
```

It defaults to SCALE, o/e, and one requested eigenvector. The solver retains two
additional internal Ritz vectors for convergence, but only `nv` vectors are
written. Eigenvalues, gaps, both conditioning ratios, and the tolerance-scaled
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
