#!/bin/bash

# Exit on error
set -e

echo "🔧 Setting up build environment..."

# Install dependencies using apt
echo "📦 Installing dependencies..."
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    gcc \
    g++ \
    libopenblas-dev \
    liblapack-dev \
    liblapacke-dev \
    libcurl4-openssl-dev \
    zlib1g-dev

# Check if straw is available
STRAW_PATH="$HOME/straw"
if [ ! -d "$STRAW_PATH" ]; then
    echo "⚠️  Warning: straw library not found at $STRAW_PATH"
    echo "Please install straw library and place it in $STRAW_PATH"
    echo "or modify this script with the correct path"
    exit 1
fi

echo "🔨 Building executables..."

# Common compiler flags
COMMON_FLAGS="-O2 -Wno-format-security -I/usr/include -I$STRAW_PATH/C++"
COMMON_LIBS="-L/usr/lib -lz -lcurl -lpthread -lopenblas -llapack -llapacke"

# First compile straw library
echo "Building straw library..."
g++ $COMMON_FLAGS -c "$STRAW_PATH/C++/straw.cpp" -o straw.o

# Compile Lan.exe
echo "Building Lan.exe..."
g++ $COMMON_FLAGS -std=c++11 -o Lan.exe \
    s_fLan.cpp \
    s_fSOLan.c \
    s_dthMul.c \
    hgFlipSign.c \
    straw.o \
    -I. \
    $COMMON_LIBS

# Compile GWev.exe
echo "Building GWev.exe..."
g++ $COMMON_FLAGS -std=c++11 -o GWev.exe \
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