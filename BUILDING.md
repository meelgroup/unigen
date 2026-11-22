# Building UniGen

Prefer the [release binaries](https://github.com/meelgroup/unigen/releases) or
Nix (`nix shell github:meelgroup/unigen#unigen`). The
[GitHub Action](https://github.com/meelgroup/unigen/actions/workflows/build.yml)
shows the exact steps for Linux, macOS and Windows.

## Dependencies

System libraries: GMP (with C++ support), MPFR, and optionally zlib (gzipped
CNFs). On Debian/Ubuntu:

```shell
sudo apt-get install libgmp-dev libmpfr-dev zlib1g-dev
```

The SAT-side dependencies are CryptoMiniSat, Arjun and ApproxMC, which bring
their own (SBVA, CaDiCaL, CaDiBack). CMake finds each one in this order:

1. a pre-built tree given as `-D<name>_DIR=...` (e.g. `-Darjun_DIR=../arjun/build`);
2. a sibling checkout next to `unigen/` (e.g. `../arjun`), built as part of UniGen;
3. otherwise it is fetched from GitHub.

So a plain build needs only the system libraries:

```shell
git clone https://github.com/meelgroup/unigen
cd unigen
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Dynamic build against pre-built dependencies

For development, keep all repositories side by side and build each one in its
`build/` directory with its `build_norm.sh`: cadical, cadiback, cryptominisat,
sbva, arjun, approxmc. Then `unigen/build/build_norm.sh` points CMake at them
(`-D<name>_DIR`), builds, and runs the tests.

## Static build

First build GMP with position-independent code. GMP's hand-optimised assembly
is normally compiled without `-fPIC` (fine for a native static binary), but
`-fPIC` is required whenever the static `.a` is linked into a shared object.
`--with-pic` makes both work:

```shell
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --enable-cxx --enable-shared --with-pic
make -j$(nproc)
sudo make install
cd ..
```

Then point CMake at the installed static GMP (use `/usr/local/lib/`, not a
custom build directory, which may be compiled for the wrong architecture):

```shell
mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=OFF \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    ..
make -j$(nproc)
```

`scripts/build_static.sh` does this against side-by-side pre-built
dependencies (which must themselves be static), and
`scripts/build_static_release.sh` also names the binary after the SHA1s of
UniGen and its dependencies.

## WebAssembly

The in-browser version is built with Emscripten. GMP and MPFR must first be
built for wasm and installed into `$EMINSTALL` (Emscripten's sysroot), and each
dependency built for wasm in its own `build-wasm/`. Then
`scripts/build_emscripten.sh` builds UniGen in `build-wasm/`, so the native
`build/` stays, and puts `unigen.js` and `unigen.wasm` into `html/`.

## Useful CMake options

| Option | Effect |
|---|---|
| `-DCMAKE_BUILD_TYPE=Debug` / `Release` / `RelWithDebInfo` | build type |
| `-DBUILD_SHARED_LIBS=OFF` | static binary |
| `-DENABLE_TESTING=ON` | build the tests (`ctest`) |
| `-DSANITIZE=ON` | Clang sanitizers |
