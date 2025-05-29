#!/bin/bash

# Exit on error
set -e

echo "🔧 Setting up build environment..."

# Assume Homebrew is installed

# Install dependencies
echo "📦 Installing dependencies..."
brew install openblas
brew install lapack
brew install gcc@13  # Latest stable GCC
brew install curl
brew install zlib

# Set up environment variables for OpenBLAS and LAPACK
export LDFLAGS="-L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/lapack/lib"
export CPPFLAGS="-I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/lapack/include"

# Check if straw is available
STRAW_PATH="../../../straw"
if [ ! -d "$STRAW_PATH" ]; then
    echo "⚠️  Warning: straw library not found at $STRAW_PATH"
    echo "Please install straw library and place it in $STRAW_PATH"
    echo "or modify this script with the correct path"
    exit 1
fi

echo "🔨 Building executables..."

# Common compiler flags
COMMON_FLAGS="-O2 -Wno-format-security -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/lapack/include -I$STRAW_PATH/C++"
COMMON_LIBS="-L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/lapack/lib -lz -lcurl -lpthread -lopenblas -llapack -llapacke"

# First compile straw library
echo "Building straw library..."
g++-13 $COMMON_FLAGS -c "$STRAW_PATH/C++/straw.cpp" -o straw.o

# Compile Lan.exe
echo "Building Lan.exe..."
g++-13 $COMMON_FLAGS -std=c++11 -o Lan.exe \
    s_fLan.cpp \
    s_fSOLan.c \
    s_dthMul.c \
    hgFlipSign.c \
    straw.o \
    -I. \
    $COMMON_LIBS

# Compile GWev.exe
echo "Building GWev.exe..."
g++-13 $COMMON_FLAGS -std=c++11 -o GWev.exe \
    s_fGW.cpp \
    getGWMatrix.cpp \
    s_fSOLan.c \
    s_dthMul.c \
    straw.o \
    $COMMON_LIBS

# Clean up object files
rm -f straw.o

echo "✅ Build completed successfully!"
echo
echo "You can now run:"
echo "  ./Lan.exe  - for chromosome-specific analysis"
echo "  ./GWev.exe - for genome-wide analysis" 