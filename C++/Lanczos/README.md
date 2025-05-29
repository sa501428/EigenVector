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

### Chromosome-Specific Analysis (Lan.exe)
```bash
./Lan.exe [options] <hicfile> <chromosome> <outbase> <resolution> [nv]

Options:
  -o           Use observed matrix instead of observed/expected (o/e) matrix
  -t <float>   Set tolerance (default: 1.0e-7)
  -e <float>   Set epsilon (default: 1.0e-8)
  -I <int>     Set maximum iterations (default: 200)
  -n <string>  Set normalization method (default: NONE)
  -T <int>     Set number of threads (default: 1)
  -v <int>     Set verbosity level (default: 1)
```

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
