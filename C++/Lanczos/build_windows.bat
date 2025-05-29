@echo off
setlocal enabledelayedexpansion

echo 🔧 Setting up build environment...

REM Check if MinGW-w64 is installed and in PATH
where g++ >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: g++ not found. Please install MinGW-w64 and add it to your PATH
    echo You can download it from: https://winlibs.com/
    exit /b 1
)

REM Set paths - MODIFY THESE AS NEEDED
set "MINGW_PATH=C:\mingw64"
set "STRAW_PATH=%USERPROFILE%\straw"
set "OPENBLAS_PATH=%MINGW_PATH%\opt\openblas"

REM Check if straw exists
if not exist "%STRAW_PATH%" (
    echo ⚠️  Warning: straw library not found at %STRAW_PATH%
    echo Please install straw library and place it in %STRAW_PATH%
    echo or modify this script with the correct path
    exit /b 1
)

echo 🔨 Building executables...

REM Common compiler flags
set "COMMON_FLAGS=-O2 -Wno-format-security -I%MINGW_PATH%\include -I%STRAW_PATH%\C++"
set "COMMON_LIBS=-L%MINGW_PATH%\lib -lz -lcurl -lpthread -lopenblas -llapack -llapacke"

REM First compile straw library
echo Building straw library...
g++ %COMMON_FLAGS% -c "%STRAW_PATH%\C++\straw.cpp" -o straw.o

REM Compile Lan.exe
echo Building Lan.exe...
g++ %COMMON_FLAGS% -std=c++11 -o Lan.exe ^
    s_fLan.cpp ^
    s_fSOLan.c ^
    s_dthMul.c ^
    hgFlipSign.c ^
    straw.o ^
    -I. ^
    %COMMON_LIBS%

REM Compile GWev.exe
echo Building GWev.exe...
g++ %COMMON_FLAGS% -std=c++11 -o GWev.exe ^
    s_fGW.cpp ^
    getGWMatrix.cpp ^
    s_fSOLan.c ^
    s_dthMul.c ^
    straw.o ^
    %COMMON_LIBS%

REM Clean up object files
del straw.o

echo ✅ Build completed successfully!
echo.
echo You can now run:
echo   Lan.exe  - for chromosome-specific analysis
echo   GWev.exe - for genome-wide analysis

endlocal 