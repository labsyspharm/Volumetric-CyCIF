#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${PARALLEL_IMARIS_WRITER_ROOT:?Set PARALLEL_IMARIS_WRITER_ROOT to an external Parallel_Imaris_Writer source checkout}"
SRC_FILE="${PARALLEL_IMARIS_WRITER_ROOT}/src/main.cpp"
OUTPUT_BIN="${1:-${REPO_ROOT}/parallelimariswriter}"

CXX_BIN="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--O3 -DNDEBUG}"
LDFLAGS_EXTRA="${LDFLAGS_EXTRA:-}"

: "${IMARIS_WRITER_ROOT:?Set IMARIS_WRITER_ROOT to the ImarisWriter install root}"
: "${CPP_TIFF_ROOT:?Set CPP_TIFF_ROOT to the cpp-tiff install root}"
: "${CPP_ZARR_ROOT:?Set CPP_ZARR_ROOT to the cpp-zarr install root}"
: "${HDF5_ROOT:?Set HDF5_ROOT to the HDF5 install root}"
: "${ZLIB_ROOT:?Set ZLIB_ROOT to the zlib install root}"

if [[ ! -f "${SRC_FILE}" ]]; then
  echo "Missing source file: ${SRC_FILE}" >&2
  exit 1
fi

if [[ ! -d "${IMARIS_WRITER_ROOT}/include" ]]; then
  echo "Missing ImarisWriter headers under ${IMARIS_WRITER_ROOT}/include" >&2
  exit 1
fi

if [[ ! -f "${IMARIS_WRITER_ROOT}/lib/libImarisWriter_static.a" ]]; then
  echo "Missing libImarisWriter_static.a under ${IMARIS_WRITER_ROOT}/lib" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUTPUT_BIN}")"

set -x
"${CXX_BIN}" \
  ${CXXFLAGS} \
  -Wl,-rpath,'$ORIGIN' \
  -I"${IMARIS_WRITER_ROOT}/include" \
  -I"${CPP_TIFF_ROOT}/include" \
  -I"${CPP_ZARR_ROOT}/include" \
  -L"${CPP_TIFF_ROOT}/lib" \
  -L"${CPP_TIFF_ROOT}/lib64" \
  -L"${CPP_ZARR_ROOT}/lib" \
  -L"${CPP_ZARR_ROOT}/lib64" \
  -lpthread \
  -fopenmp \
  -static-libstdc++ \
  -static-libgcc \
  -ldl \
  -lcppTiff \
  -lcppZarr \
  "${SRC_FILE}" \
  -o "${OUTPUT_BIN}" \
  "${IMARIS_WRITER_ROOT}/lib/libImarisWriter_static.a" \
  "${HDF5_ROOT}/lib/libhdf5.a" \
  "${ZLIB_ROOT}/lib/libz.a" \
  ${LDFLAGS_EXTRA}
set +x

echo
echo "Built: ${OUTPUT_BIN}"
echo "Next: run ldd ${OUTPUT_BIN} and then test on one dataset."
