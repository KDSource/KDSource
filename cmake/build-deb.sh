#!/bin/bash
# Build the native (non-Python) KDSource .deb with CMake + CPack.
# Needs: build-essential cmake libxml2-dev dpkg-dev, and MCPL (findable via
# mcpl-config or -DMCPL_DIR=...).
# Usage: cmake/build-deb.sh [builddir] [extra cmake args...]   (default builddir: <repo>/tmp)
# e.g.   cmake/build-deb.sh tmp -DCPACK_DEBIAN_PACKAGE_DEPENDS="libxml2, libmcpl (>= 2.2.8)"
set -euo pipefail
src="$(cd "$(dirname "$0")/.." && pwd)"
build="${1:-$src/tmp}"
shift || true
cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DKDS_ENABLE_CPACK=ON "$@"
cmake --build "$build" -j"$(nproc)"
(cd "$build" && cpack -G DEB)
dpkg-deb -I "$build"/kdsource_*.deb
