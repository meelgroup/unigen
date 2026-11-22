[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![build](https://github.com/meelgroup/unigen/workflows/build/badge.svg)

# UniGen: Almost-Uniform Sampler
UniGen is the state of the art almost-uniform sampler  utilizing an improved
version of CryptoMiniSat to handle problems of size and complexity that were
not possible before. The current version is based on work Mate Soos, Stephan
Gocht, and Kuldeep S. Meel, as [published in
CAV-20](http://www.cs.toronto.edu/~meel/Papers/cav20-sgm.pdf). Please see below
for credits.  A large part of the work is in CryptoMiniSat
[here](https://github.com/msoos/cryptominisat).

## Installation
You can try out UniGen [from your browser](https://www.msoos.org/unigen/).

It is strongly recommended to not build, but to use the precompiled
binaries as in our [release](https://github.com/meelgroup/unigen/releases).
The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:
```shell
nix shell github:meelgroup/unigen#unigen
```

Then you will have `unigen` binary available and ready to use.

### Building from source
A plain CMake build fetches the dependencies (CryptoMiniSat, Arjun, ApproxMC
and theirs) from GitHub, or uses them from sibling directories if they are
checked out next to UniGen:
```shell
mkdir build && cd build
cmake ..
make -j$(nproc)
```

With the dependencies already built, `scripts/build_norm.sh` (run from
`build/`) builds against their `build/` directories. To build a static binary
(`scripts/build_static.sh`), you first need GMP's static libraries:
```shell
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --enable-cxx --enable-shared --with-pic
make -j8
sudo make install
cd ..
```

## Sampling/Projection Set
For some applications, one is not interested in solutions over all the
variables and instead interested in sampling over the solutions projected to a
subset of variables, called sampling set. UniGen allows you to specify the
sampling set using the following modified version of DIMACS format:
```plain
$ cat myfile.cnf
c p show 1 3 4 6 7 8 10 0
p cnf 500 1
3 4 0
```
Above, using the `c p show` line, we declare that only variables 1, 3, 4, 6, 7,
8 and 10 form part of the sampling set out of the CNF's 500 variables
`1,2...500`. This line must end with a 0. The legacy directive `c ind ... 0`
declares the same thing.  Naturally, if your sampling set only
contains 7 variables, then the maximum number of solutions restricted to
sampling set can only be at most 2^7 = 128. This is true even if your CNF has
thousands of variables.

UniGen samples over the CNF that [Arjun](https://github.com/meelgroup/arjun)
simplified, and maps each sample back to the input with Arjun's definitions of
the variables it removed. By default only the sampling set is printed; with
`--fullsol 1` every variable of the input gets a value, i.e. each sample is a
full solution of the input CNF.

## Library Use
Below is an example library use:
```cpp
#include <unigen/unigen.h>
#include <approxmc/approxmc.h>
#include <arjun/arjun.h>
#include <iostream>
#include <memory>
#include <vector>
using std::vector;
using namespace CMSat;

void mycallback(const vector<int>& solution, void*) {
    for (const int lit: solution) std::cout << lit << " ";
    std::cout << "0" << std::endl;
}

int main() {
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpz>();
    ApproxMC::AppMC appmc(fg);
    UniGen::UniG unigen(&appmc);
    unigen.set_callback(mycallback, nullptr);

    appmc.new_vars(3);
    appmc.add_clause({Lit(0, true), Lit(1, true)});
    appmc.set_sampl_vars({0, 1, 2});

    auto sol_count = appmc.count();
    unigen.sample(&sol_count, 5);
    return 0;
}
```

## Guarantees
UniGen ensures that the generated distribution is within
(1+\epsilon)-multiplicative factor of the ideal uniform distribution.

## Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new
issue](https://github.com/meelgroup/unigen/issues/new). All issues are
responded to promptly.

## How to Cite
If you use UniGen, please cite the following papers:
[CAV20](https://www.cs.toronto.edu/~meel/bib/SGM20.bib),
[TACAS15](https://www.cs.toronto.edu/~meel/bib/CFMSV15a.bib), and
[DAC14](https://www.cs.toronto.edu/~meel/bib/CMV14.bib).

UniGen builds on a series of papers on hashing-based approach: [Related
Publications](https://www.cs.toronto.edu/~meel/publications.html)

## Contributors
UniGen is based on research papers published over the period of 2012-20 and
co-authored by (in alphabetical order): Supratik Chakraborty, Daniel Fremont,
Stephan Gocht, Kuldeep S. Meel, Sanjit A. Seshia, Mate Soos, and Moshe Y.
Vardi.
