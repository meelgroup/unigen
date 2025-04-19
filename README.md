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

## Sampling Set
For some applications, one is not interested in solutions over all the
variables and instead interested in sampling over the solutions projected to a
subset of variables, called sampling set. UniGen allows you to specify the
sampling set using the following modified version of DIMACS format:
```plain
$ cat myfile.cnf
c ind 1 3 4 6 7 8 10 0
p cnf 500 1
3 4 0
```
Above, using the `c ind` line, we declare that only variables 1, 3, 4, 6, 7, 8
and 10 form part of the sampling set out of the CNF's 500 variables
`1,2...500`. This line must end with a 0.  Naturally, if your sampling set only
contains 7 variables, then the maximum number of solutions restricted to
sampling set can only be at most 2^7 = 128. This is true even if your CNF has
thousands of variables.

## Building
It is highly discouraged to build UniGen from source. Instead, please use the
provided static binaries for your platform from our
[releases](https://github.com/meelgroup/unigen/releases/tag/2.5.7). In case you
must compile from source, please follow the [GitHub
action](https://github.com/meelgroup/unigen/actions/workflows/build.yml), whih
builds everything from source for all supported platforms.

## Library Use
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
    lits.push_back(Lit(2, true));https://github.com/meelgroup/unigen
    appmc->addClause(lits);

    sol_count = appmc->count();
    unigen->sample(&sol_count, num_samples);
    return 0;
}
```

### Guarantees
UniGen ensures that the generated distribution is within
(1+\epsilon)-multiplicative factor of the ideal uniform distribution.


### Issues, questions, bugs, etc.
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
