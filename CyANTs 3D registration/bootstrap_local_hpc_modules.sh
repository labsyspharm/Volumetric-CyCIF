#!/usr/bin/env bash
# Developed by Alex Wong
# Cite: Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in
# Preprint: https://doi.org/10.64898/2026.05.17.725158
# Registration method: ANTsX/ANTsPy - https://github.com/ANTsX/ANTsPy

set -euo pipefail

PREFIX="${1:-$HOME/opt/hpc-bootstrap}"
PREFIX="$(mkdir -p "$PREFIX" && cd "$PREFIX" && pwd)"
SOFTWARE_ROOT="${PREFIX}/software"
MODULE_ROOT="${PREFIX}/modulefiles"
SRC_ROOT="${PREFIX}/src"
BUILD_ROOT="${PREFIX}/build"

ZLIB_VERSION="${ZLIB_VERSION:-1.3.2}"
LZ4_VERSION="${LZ4_VERSION:-1.10.0}"
HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
CMAKE_VERSION="${CMAKE_VERSION:-4.2.3}"

ZLIB_URL="${ZLIB_URL:-https://zlib.net/current/zlib.tar.gz}"
LZ4_URL="${LZ4_URL:-https://github.com/lz4/lz4/archive/refs/tags/v${LZ4_VERSION}.tar.gz}"
HDF5_URL="${HDF5_URL:-https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_${HDF5_VERSION}.tar.gz}"
CMAKE_URL="${CMAKE_URL:-https://cmake.org/files/v4.2/cmake-${CMAKE_VERSION}.tar.gz}"

BUILD_JOBS="${BUILD_JOBS:-8}"

CC_BIN="${CC:-$(command -v cc || true)}"
CXX_BIN="${CXX:-$(command -v c++ || true)}"
if [[ -z "${CC_BIN}" || -z "${CXX_BIN}" ]]; then
  echo "Need a bootstrap C/C++ compiler (cc and c++)." >&2
  exit 1
fi

DOWNLOAD_TOOL=""
if command -v curl >/dev/null 2>&1; then
  DOWNLOAD_TOOL="curl"
elif command -v wget >/dev/null 2>&1; then
  DOWNLOAD_TOOL="wget"
else
  echo "Need curl or wget available to download sources." >&2
  exit 1
fi

mkdir -p "${SOFTWARE_ROOT}" "${MODULE_ROOT}" "${SRC_ROOT}" "${BUILD_ROOT}"

fetch() {
  local url="$1"
  local dest="$2"
  if [[ -f "${dest}" ]]; then
    return 0
  fi
  echo "[download] ${url}"
  if [[ "${DOWNLOAD_TOOL}" == "curl" ]]; then
    curl -L --fail --output "${dest}" "${url}"
  else
    wget -O "${dest}" "${url}"
  fi
}

extract_tarball() {
  local tarball="$1"
  local dest_root="$2"
  local marker="$3"
  local dest_dir="${dest_root}/${marker}"
  rm -rf "${dest_dir}"
  mkdir -p "${dest_dir}"
  tar -xf "${tarball}" -C "${dest_dir}" --strip-components=1
  printf '%s\n' "${dest_dir}"
}

prepend_path() {
  local var_name="$1"
  local value="$2"
  local current="${!var_name:-}"
  if [[ -n "${current}" ]]; then
    printf '%s:%s' "${value}" "${current}"
  else
    printf '%s' "${value}"
  fi
}

write_modulefile() {
  local name="$1"
  local version="$2"
  local prefix="$3"
  local extra_body="${4:-}"
  local moddir="${MODULE_ROOT}/${name}"
  mkdir -p "${moddir}"
  cat > "${moddir}/${version}" <<EOF
#%Module1.0
proc ModulesHelp { } {
    puts stderr "${name} ${version}"
}
module-whatis "${name} ${version}"
set root ${prefix}
prepend-path PATH \$root/bin
prepend-path CPATH \$root/include
prepend-path LIBRARY_PATH \$root/lib
prepend-path LIBRARY_PATH \$root/lib64
prepend-path LD_LIBRARY_PATH \$root/lib
prepend-path LD_LIBRARY_PATH \$root/lib64
prepend-path PKG_CONFIG_PATH \$root/lib/pkgconfig
prepend-path PKG_CONFIG_PATH \$root/lib64/pkgconfig
setenv ${name^^}_ROOT \$root
${extra_body}
EOF
}

write_system_gcc_modulefile() {
  local gcc_bin
  gcc_bin="$(command -v gcc || true)"
  if [[ -z "${gcc_bin}" ]]; then
    echo "[warn] gcc not found; skipping gcc modulefile."
    return 0
  fi
  local gcc_root
  gcc_root="$(cd "$(dirname "${gcc_bin}")/.." && pwd)"
  local moddir="${MODULE_ROOT}/gcc"
  mkdir -p "${moddir}"
  cat > "${moddir}/system" <<EOF
#%Module1.0
proc ModulesHelp { } {
    puts stderr "System GCC compiler"
}
module-whatis "System GCC compiler"
set root ${gcc_root}
prepend-path PATH \$root/bin
prepend-path LIBRARY_PATH \$root/lib
prepend-path LIBRARY_PATH \$root/lib64
prepend-path LD_LIBRARY_PATH \$root/lib
prepend-path LD_LIBRARY_PATH \$root/lib64
setenv GCC_ROOT \$root
EOF
}

build_zlib() {
  local install_prefix="${SOFTWARE_ROOT}/zlib/${ZLIB_VERSION}"
  if [[ -x "${install_prefix}/bin/zlib-flate" || -f "${install_prefix}/lib/libz.a" || -f "${install_prefix}/lib/libz.so" ]]; then
    echo "[skip] zlib ${ZLIB_VERSION}"
    write_modulefile "zlib" "${ZLIB_VERSION}" "${install_prefix}"
    return 0
  fi
  local tarball="${SRC_ROOT}/zlib-${ZLIB_VERSION}.tar.gz"
  fetch "${ZLIB_URL}" "${tarball}"
  local src_dir
  src_dir="$(extract_tarball "${tarball}" "${BUILD_ROOT}" "zlib-${ZLIB_VERSION}")"
  (
    cd "${src_dir}"
    ./configure --prefix="${install_prefix}"
    make -j"${BUILD_JOBS}"
    make install
  )
  write_modulefile "zlib" "${ZLIB_VERSION}" "${install_prefix}"
}

build_lz4() {
  local install_prefix="${SOFTWARE_ROOT}/lz4/${LZ4_VERSION}"
  if [[ -x "${install_prefix}/bin/lz4" || -f "${install_prefix}/lib/liblz4.a" || -f "${install_prefix}/lib/liblz4.so" ]]; then
    echo "[skip] lz4 ${LZ4_VERSION}"
    write_modulefile "lz4" "${LZ4_VERSION}" "${install_prefix}"
    return 0
  fi
  local tarball="${SRC_ROOT}/lz4-${LZ4_VERSION}.tar.gz"
  fetch "${LZ4_URL}" "${tarball}"
  local src_dir
  src_dir="$(extract_tarball "${tarball}" "${BUILD_ROOT}" "lz4-${LZ4_VERSION}")"
  (
    cd "${src_dir}"
    make -j"${BUILD_JOBS}"
    make install PREFIX="${install_prefix}"
  )
  write_modulefile "lz4" "${LZ4_VERSION}" "${install_prefix}"
}

build_cmake() {
  local install_prefix="${SOFTWARE_ROOT}/cmake/${CMAKE_VERSION}"
  if [[ -x "${install_prefix}/bin/cmake" ]]; then
    echo "[skip] cmake ${CMAKE_VERSION}"
    write_modulefile "cmake" "${CMAKE_VERSION}" "${install_prefix}"
    return 0
  fi
  local tarball="${SRC_ROOT}/cmake-${CMAKE_VERSION}.tar.gz"
  fetch "${CMAKE_URL}" "${tarball}"
  local src_dir
  src_dir="$(extract_tarball "${tarball}" "${BUILD_ROOT}" "cmake-${CMAKE_VERSION}")"
  (
    cd "${src_dir}"
    ./bootstrap --prefix="${install_prefix}" CC="${CC_BIN}" CXX="${CXX_BIN}"
    make -j"${BUILD_JOBS}"
    make install
  )
  write_modulefile "cmake" "${CMAKE_VERSION}" "${install_prefix}"
}

build_hdf5() {
  local install_prefix="${SOFTWARE_ROOT}/hdf5/${HDF5_VERSION}"
  if [[ -x "${install_prefix}/bin/h5dump" || -f "${install_prefix}/lib/libhdf5.a" || -f "${install_prefix}/lib/libhdf5.so" ]]; then
    echo "[skip] hdf5 ${HDF5_VERSION}"
    write_modulefile "hdf5" "${HDF5_VERSION}" "${install_prefix}"
    return 0
  fi
  local tarball="${SRC_ROOT}/hdf5-${HDF5_VERSION}.tar.gz"
  fetch "${HDF5_URL}" "${tarball}"
  local src_dir
  src_dir="$(extract_tarball "${tarball}" "${BUILD_ROOT}" "hdf5-${HDF5_VERSION}")"
  local zlib_prefix="${SOFTWARE_ROOT}/zlib/${ZLIB_VERSION}"
  (
    cd "${src_dir}"
    export CPPFLAGS="-I${zlib_prefix}/include"
    export LDFLAGS="-L${zlib_prefix}/lib -L${zlib_prefix}/lib64"
    ./configure \
      --prefix="${install_prefix}" \
      --with-zlib="${zlib_prefix}" \
      --enable-build-mode=production \
      --disable-dependency-tracking
    make -j"${BUILD_JOBS}"
    make install
  )
  write_modulefile "hdf5" "${HDF5_VERSION}" "${install_prefix}"
}

cat <<EOF
[info] Installing user-space HPC software into:
  ${PREFIX}
[info] Bootstrap compiler:
  CC=${CC_BIN}
  CXX=${CXX_BIN}
[info] Planned versions:
  zlib=${ZLIB_VERSION}
  lz4=${LZ4_VERSION}
  hdf5=${HDF5_VERSION}
  cmake=${CMAKE_VERSION}
EOF

build_zlib
build_lz4
build_cmake
build_hdf5
write_system_gcc_modulefile

cat <<EOF

[done] Installed software under:
  ${SOFTWARE_ROOT}

[done] Modulefiles written under:
  ${MODULE_ROOT}

Next shell steps:
  module use ${MODULE_ROOT}
  module load gcc/system
  module load zlib/${ZLIB_VERSION}
  module load lz4/${LZ4_VERSION}
  module load cmake/${CMAKE_VERSION}
  module load hdf5/${HDF5_VERSION}

Then you can run:
  conda activate cyants
  export HDF5_ROOT=${SOFTWARE_ROOT}/hdf5/${HDF5_VERSION}
  export ZLIB_ROOT=${SOFTWARE_ROOT}/zlib/${ZLIB_VERSION}
  export LZ4_ROOT=${SOFTWARE_ROOT}/lz4/${LZ4_VERSION}
  bash /Users/alexwong/Documents/Codex/cyants/setup_cyants_imaris.sh \$HOME/src/ImarisWriter
EOF
