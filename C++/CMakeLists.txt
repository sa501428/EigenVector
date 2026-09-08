cmake_minimum_required(VERSION 3.16)
project(possum_solvers LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(STRAW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../hic-infrastructure/straw/C++"
    CACHE PATH "Directory containing current straw C++ sources")
if(NOT EXISTS "${STRAW_DIR}/straw.h" OR NOT EXISTS "${STRAW_DIR}/straw.cpp")
    message(FATAL_ERROR "STRAW_DIR must contain straw.h and straw.cpp; got ${STRAW_DIR}")
endif()

find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_path(ZSTD_INCLUDE_DIR zstd.h
    HINTS /opt/homebrew/opt/zstd/include /usr/local/include /usr/include)
find_library(ZSTD_LIBRARY NAMES zstd
    HINTS /opt/homebrew/opt/zstd/lib /usr/local/lib /usr/lib)
find_path(LAPACKE_INCLUDE_DIR lapacke.h
    HINTS /opt/homebrew/opt/lapack/include /usr/local/include /usr/include)
find_library(LAPACKE_LIBRARY NAMES lapacke
    HINTS /opt/homebrew/opt/lapack/lib /usr/local/lib /usr/lib)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
    message(FATAL_ERROR "zstd headers/library not found")
endif()
if(NOT LAPACKE_INCLUDE_DIR OR NOT LAPACKE_LIBRARY)
    message(FATAL_ERROR "LAPACKE headers/library not found")
endif()

add_library(possum_straw STATIC
    "${STRAW_DIR}/straw.cpp"
    "${STRAW_DIR}/straw_v10.cpp"
    "${STRAW_DIR}/hbs.cpp")
target_include_directories(possum_straw PUBLIC "${STRAW_DIR}" PRIVATE "${ZSTD_INCLUDE_DIR}")
target_link_libraries(possum_straw PUBLIC CURL::libcurl ZLIB::ZLIB "${ZSTD_LIBRARY}" Threads::Threads)
target_compile_options(possum_straw PRIVATE -w)

set(POWER_COMMON
    PowerMethod/theBestEigen.c PowerMethod/thdMul.c PowerMethod/hgFlipSign.c)
set_source_files_properties(${POWER_COMMON} PowerMethod/mainEigen.c PROPERTIES LANGUAGE CXX)

add_executable(possum_power_chr PowerMethod/theEigenVector_flip_new.cpp ${POWER_COMMON})
set_target_properties(possum_power_chr PROPERTIES OUTPUT_NAME ev.exe)
target_include_directories(possum_power_chr PRIVATE PowerMethod)
target_link_libraries(possum_power_chr PRIVATE possum_straw Threads::Threads)

add_executable(possum_power_genome PowerMethod/GWevIntra_new.cpp ${POWER_COMMON})
set_target_properties(possum_power_genome PROPERTIES OUTPUT_NAME GWevIntra.exe)
target_include_directories(possum_power_genome PRIVATE PowerMethod)
target_link_libraries(possum_power_genome PRIVATE possum_straw Threads::Threads)

add_executable(possum_power_ascii
    PowerMethod/mainEigen.c PowerMethod/theBestEigen.c PowerMethod/thdMul.c)
set_target_properties(possum_power_ascii PROPERTIES OUTPUT_NAME mainEigen.exe)
target_include_directories(possum_power_ascii PRIVATE PowerMethod)
target_link_libraries(possum_power_ascii PRIVATE Threads::Threads)

set(LANCZOS_COMMON
    Lanczos/s_fSOLan.c Lanczos/s_dthMul.c Lanczos/hgFlipSign.c)
set_source_files_properties(${LANCZOS_COMMON} PROPERTIES LANGUAGE CXX)

add_executable(possum2_chr Lanczos/s_fLan.cpp ${LANCZOS_COMMON})
set_target_properties(possum2_chr PROPERTIES OUTPUT_NAME LanChr.exe)
target_include_directories(possum2_chr PRIVATE Lanczos "${LAPACKE_INCLUDE_DIR}")
target_link_libraries(possum2_chr PRIVATE possum_straw "${LAPACKE_LIBRARY}" ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(possum2_genome Lanczos/LanGenome.cpp)
set_target_properties(possum2_genome PROPERTIES OUTPUT_NAME Lan.exe)
target_include_directories(possum2_genome PRIVATE Lanczos)
target_link_libraries(possum2_genome PRIVATE possum_straw)

add_executable(possum2_trans
    Lanczos/s_fGW.cpp Lanczos/getGWMatrix.cpp ${LANCZOS_COMMON})
set_target_properties(possum2_trans PROPERTIES OUTPUT_NAME GWev.exe)
target_include_directories(possum2_trans PRIVATE Lanczos "${LAPACKE_INCLUDE_DIR}")
target_link_libraries(possum2_trans PRIVATE possum_straw "${LAPACKE_LIBRARY}" ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
