#!/usr/bin/env bash
# Builds unigen for wasm in its own build-wasm/ (the native build/ stays) and
# installs it into $EMINSTALL. Run from anywhere, after source ~/emscripten.sh.
set -euo pipefail
REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
SAT_DIR="$(cd "$REPO/.." && pwd)"
mkdir -p "$REPO/build-wasm"
cd "$REPO/build-wasm"
emcmake cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX="${EMINSTALL}" \
    -Dcadical_DIR="${SAT_DIR}/cadical/build-wasm" \
    -Dcadiback_DIR="${SAT_DIR}/cadiback/build-wasm" \
    -Dcryptominisat5_DIR="${SAT_DIR}/cryptominisat/build-wasm" \
    -Dsbva_DIR="${SAT_DIR}/sbva/build-wasm" \
    -Darjun_DIR="${SAT_DIR}/arjun/build-wasm" \
    -Dapproxmc_DIR="${SAT_DIR}/approxmc/build-wasm" \
    "$REPO"
emmake make -j"$(nproc)"
emmake make install
cp unigen.wasm "$REPO/html"
cp "${EMINSTALL}/bin/unigen.js" "$REPO/html"
