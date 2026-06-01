#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 /absolute/path/to/ImarisWriter [build_dir]" >&2
  exit 1
fi

IMARISWRITER_SRC="$(cd "$1" && pwd)"
BUILD_DIR="${2:-${IMARISWRITER_SRC}/release}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "Activate your cyants conda env first, for example: conda activate cyants" >&2
  exit 1
fi

if [[ ! -f "${IMARISWRITER_SRC}/CMakeLists.txt" ]]; then
  echo "Could not find ImarisWriter source at ${IMARISWRITER_SRC}" >&2
  exit 1
fi

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake \
  ${HDF5_ROOT:+-DHDF5_ROOT:PATH=${HDF5_ROOT}} \
  ${ZLIB_ROOT:+-DZLIB_ROOT:PATH=${ZLIB_ROOT}} \
  ${LZ4_ROOT:+-DLZ4_ROOT:PATH=${LZ4_ROOT}} \
  ..

make -j"${BUILD_JOBS:-8}"

ACTIVATE_DIR="${CONDA_PREFIX}/etc/conda/activate.d"
DEACTIVATE_DIR="${CONDA_PREFIX}/etc/conda/deactivate.d"
mkdir -p "${ACTIVATE_DIR}" "${DEACTIVATE_DIR}"

ACTIVATE_FILE="${ACTIVATE_DIR}/imariswriter.sh"
DEACTIVATE_FILE="${DEACTIVATE_DIR}/imariswriter.sh"

cat > "${ACTIVATE_FILE}" <<EOF
export IMARISWRITER_ROOT="${IMARISWRITER_SRC}"
export PYTHONPATH="\${IMARISWRITER_ROOT}/python:\${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${BUILD_DIR}:\${LD_LIBRARY_PATH:-}"
EOF

cat > "${DEACTIVATE_FILE}" <<EOF
if [[ "\${PYTHONPATH:-}" == "${IMARISWRITER_SRC}/python:"* ]]; then
  export PYTHONPATH="\${PYTHONPATH#"${IMARISWRITER_SRC}/python:"}"
fi
if [[ "\${LD_LIBRARY_PATH:-}" == "${BUILD_DIR}:"* ]]; then
  export LD_LIBRARY_PATH="\${LD_LIBRARY_PATH#"${BUILD_DIR}:"}"
fi
unset IMARISWRITER_ROOT
EOF

echo
echo "ImarisWriter build finished."
echo "Reactivate the env:"
echo "  conda deactivate && conda activate $(basename "${CONDA_PREFIX}")"
echo
echo "Quick test:"
echo "  python -c \"import PyImarisWriter; print('PyImarisWriter import ok')\""
