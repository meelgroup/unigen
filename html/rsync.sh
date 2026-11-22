#!/bin/bash
set -euxo pipefail

rsync -vaP unigen.js unigen.wasm  index.html msoos.org:/var/www/unigen/
