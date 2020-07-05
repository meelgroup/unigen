#!/bin/bash

set -e

rm -rf lib* Test* tests* include tests scalmc* CM* cmake* approxmc* src-u* unig*
cmake -DENABLE_TESTING=ON ..
make -j6
make test
