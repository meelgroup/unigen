#!/usr/bin/env bash

set -euo pipefail

SOLVERS_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
echo "solvers dir: $SOLVERS_DIR"

rm -rf .cmake
rm -rf lib*
rm -rf Test*
rm -rf tests*
rm -rf include
rm -rf tests
rm -rf CM*
rm -rf cmake*
rm -rf unigen
rm -rf src-unigen
rm -rf Makefile
rm -rf deps
rm -rf _deps
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSTATIC_BINARY=ON \
    -DGMP_LIBRARY=/usr/local/lib/libgmp.a \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    -Dcadical_DIR="${SOLVERS_DIR}/cadical/build" \
    -Dcadiback_DIR="${SOLVERS_DIR}/cadiback/build" \
    -Dcryptominisat5_DIR="${SOLVERS_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SOLVERS_DIR}/sbva/build" \
    -Darjun_DIR="${SOLVERS_DIR}/arjun/build" \
    -Dapproxmc_DIR="${SOLVERS_DIR}/approxmc/build" \
    ..
make -j$(nproc)
strip unigen

HASHES=$(./unigen --version 2>&1 || true)
sha() { echo "$HASHES" | grep "c o $1 SHA1:" | grep -oP '[0-9a-f]{40}' | cut -c1-8; }
DEST="unigen_$(sha UniGen)_$(sha ApproxMC)_$(sha Arjun)_$(sha CMS)_$(sha CadiBack)"
cp unigen "$DEST"
echo "Copied to $DEST"
