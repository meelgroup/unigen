#!/bin/bash

set -e
rm -rf cm* CM* lib* unigen* Testing* tests* include tests utils Make*
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL -DENABLE_TESTING=OFF ..
emmake make -j26
emmake make install
cp unigen.wasm ../html
cp $EMINSTALL/bin/unigen.js ../html
