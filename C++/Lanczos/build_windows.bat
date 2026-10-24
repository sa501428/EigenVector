@echo off
setlocal
REM Run from a MinGW-w64 environment with CMake, curl, zlib, zstd,
REM pthreads, BLAS, LAPACK and LAPACKE installed.
REM Set STRAW_DIR to the directory containing straw.h and straw.cpp.
set "SOURCE=%~dp0.."
set "BUILD=%~dp0..\build"
set "STRAW_OPTION="
if defined STRAW_DIR set STRAW_OPTION=-DSTRAW_DIR="%STRAW_DIR%"
cmake -S "%SOURCE%" -B "%BUILD%" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="%~dp0." %STRAW_OPTION% %*
if errorlevel 1 exit /b 1
cmake --build "%BUILD%" --parallel
exit /b %errorlevel%
