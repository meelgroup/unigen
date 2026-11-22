# pyunigen: bindings to the UniGen almost uniform sampler

This directory provides Python bindings to UniGen on the C++ level,
i.e. when importing pyunigen, UniGen becomes part of the
Python process itself.

## Installing

```
pip install pyunigen
```

Pre-built wheels are available for Linux (x86-64, ARM64) and macOS (Apple
Silicon).

## Compiling
If you don't want to use the pip package, you can build it from a checkout.
It needs GMP and MPFR; CMake builds UniGen's other dependencies from sibling
checkouts next to `unigen/`, or fetches them from GitHub:

```
pip install .
```

## Usage

The `pyunigen` module has one object, `Sampler` that has two functions
`sample` and `add_clause`.

The function `add_clause()` takes an iterable list of literals such as
`[1, 2]` which represents the truth `1 or 2 = True`. For example,
`add_clause([1])` sets variable `1` to `True`.

The function `sample()` samples the system of equations that have been added
with `add_clause()`:

```
>>> from pyunigen import Sampler
>>> c = Sampler()
>>> c.add_clause([1, 5])
>>> c.add_clause([10, 11, 12])
>>> cells, hashes, samples = c.sample(num=2, sampling_set=range(1,5))
>>> print("There are approx. ", cells*2**hashes, " solutions over the sampling set. Samples: ", samples)
There are approx.  16  solutions over the sampling set. Samples:  [[1, -2, 3, -4], [1, 2, -3, -4]]
```

The return value is a tuple of cells and hashes. Which gives how many solutions
there are, probabilistically approximately

The samples are over the sampling set, by default all variables. UniGen
simplifies the CNF with [Arjun](https://github.com/meelgroup/arjun) and maps
every sample back to its variables. On an unsatisfiable CNF, `sample()` returns
`(0, 0, [])`. `sample()` can only be called once per `Sampler`.

You can give the following arguments to `Sampler`:
* `seed` -- sets the random seed
* `verbosity` -- sets the verbosity of the system (default = 0)
* `epsilon` -- Tolerance parameter of the count, i.e. sets how approximate it is. Default = 0.8
* `delta` -- Confidence parameter of the count, i.e. sets how probabilistically correct it is. Default = 0.20
* `kappa` -- Uniformity parameter (see the TACAS-15 paper). Default = 0.638
* `multisample` -- Return multiple samples from each cell. Default = False
* `arjun` -- Simplify the CNF with Arjun before sampling. Default = True
