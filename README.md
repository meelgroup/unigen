[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![build](https://github.com/meelgroup/unigen/workflows/build/badge.svg)
[![Docker Hub](https://img.shields.io/badge/docker-latest-blue.svg)](https://hub.docker.com/r/msoos/unigen/)

# UniGen3: Almost-Uniform Sampler
UniGen3 is the state of the art almost-uniform sampler  utilizing an improved version of CryptoMiniSat to handle problems of size and complexity that were not possible before. The current version is based on work Mate Soos, Stephan Gocht, and Kuldeep S. Meel, as [published in CAV-20](http://comp.nus.edu.sg/~meel/Papers/cav20-sgm.pdf). Please see below for credits.  A large part of the work is in CryptoMiniSat [here](https://github.com/msoos/cryptominisat).

## How to use the Python package

```
pip install pyunigen
python
>>> from pyunigen import Sampler
>>> c = Sampler()
>>> c.add_clause([1, 5])
>>> c.add_clause([10, 11, 12])
>>> cells, hashes, samples = c.sample(num=2, sampling_set=range(1,5))
>>> print("There are approx. ", cells*2**hashes, " solutions over the sampling set. Samples: ", samples)
There are approx. 16  solutions over the sampling set. Samples:  [[1, -2, 3, -4], [1, 2, -3, -4]]
```

## How to Build and Use on Linux:

```
sudo apt-get install build-essential cmake
sudo apt-get install zlib1g-dev libboost-program-options-dev libboost-serialization-dev
git clone https://github.com/msoos/cryptominisat
cd cryptominisat
mkdir build && cd build
cmake -DUSE_GAUSS=ON ..
make
sudo make install
cd ../..
git clone https://github.com/meelgroup/approxmc/
cd approxmc
mkdir build && cd build
cmake ..
make
sudo make install
cd ../..
git clone https://github.com/meelgroup/unigen/
cd unigen
mkdir build && cd build
cmake ..
make
sudo make install
```

If you need a static binary, you must build CryptoMiniSat, ApproxMC, and Unigen **all** with the `cmake` argument `-DSTATICCOMPILE=ON`.

### Sampling Set

For some applications, one is not interested in solutions over all the variables and instead interested in sampling over the solutions projected to a subset of variables, called sampling set. UniGen allows you to specify the sampling set using the following modified version of DIMACS format:

```
$ cat myfile.cnf
c ind 1 3 4 6 7 8 10 0
p cnf 500 1
3 4 0
```
Above, using the `c ind` line, we declare that only variables 1, 3, 4, 6, 7, 8 and 10 form part of the sampling set out of the CNF's 500 variables `1,2...500`. This line must end with a 0.  Naturally, if your sampling set only contains 7 variables, then the maximum number of solutions restricted to sampling set can only be at most 2^7 = 128. This is true even if your CNF has thousands of variables.

### Independent set
For most applications, we are want all solutions to the problem. To do this, you need to use the [MIS](https://github.com/meelgroup/mis) tool to find a small independent set of variables to your CNF. For example, for `formula.cnf` we can do:

```
docker run --rm -v `pwd`/formula.cnf:/in msoos/mis --timeout 300 /in
[...]
** Copy-paste the following line in the top of your CNF for UniGen **
c ind 3 4 7 8 10 11 14 17 18 26 30 35 36 39 42 47 60 62 67 0
```

You must copy the line starting with `c ind ...` to the top of your CNF before running ApproxMC.

### Library Use
Below is an example library use:

```
#include <unigen/unigen.h>
using namespace CMSat;
using namespace UniGen;
using std::cout;
using std::endl;

void mycallback(const std::vector<int>& solution, void*)
{
    for(uint32_t i = 0; i < solution.size(); i++) {
        cout << solution[i] <<  " ";
    }
     cout << "0" << endl;
}

int main () {
    auto appmc = new ApproxMC::AppMC;
    auto unigen = new UniG(appmc);
    appmc->set_verbosity(verbosity);
    unigen->set_callback(mycallback, NULL);
    vector<Lit> lits;

    appmc->add_variables(3);
    lits.clear();
    lits.push_back(Lit(0, true));
    lits.push_back(Lit(1, true));
    appmc->addClause(lits);

    auto sol_count = appmc->count();
    unigen->sample(&sol_count, num_samples);

    lits.clear();
    lits.push_back(Lit(0, true));
    lits.push_back(Lit(2, true));
    appmc->addClause(lits);

    sol_count = appmc->count();
    unigen->sample(&sol_count, num_samples);
    return 0;
}
```


### Guarantees
UniGen ensures that the generated distribution is within (1+\epsilon)-multiplicative factor of the ideal uniform distribution. 


### Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/unigen/issues/new). All issues are responded to promptly.

## How to Cite
If you use UniGen, please cite the following papers: [CAV20](https://www.comp.nus.edu.sg/~meel/bib/SGM20.bib), [TACAS15](https://www.comp.nus.edu.sg/~meel/bib/CFMSV15a.bib), and [DAC14](https://www.comp.nus.edu.sg/~meel/bib/CMV14.bib).

UniGen builds on a series of papers on hashing-based approach: [Related Publications](https://www.comp.nus.edu.sg/~meel/publications.html)

## Contributors
UniGen is based on research papers published over the period of 2012-20 and co-authored by (in alphabetical order): Supratik Chakraborty, Daniel Fremont, Stephan Gocht, Kuldeep S. Meel, Sanjit A. Seshia, Mate Soos, and Moshe Y. Vardi. 


