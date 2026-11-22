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
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON \
    -Dcadical_DIR="${SOLVERS_DIR}/cadical/build" \
    -Dcadiback_DIR="${SOLVERS_DIR}/cadiback/build" \
    -Dcryptominisat5_DIR="${SOLVERS_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SOLVERS_DIR}/sbva/build" \
    -Darjun_DIR="${SOLVERS_DIR}/arjun/build" \
    -Dapproxmc_DIR="${SOLVERS_DIR}/approxmc/build" \
    ..
make -j$(nproc)
make test
