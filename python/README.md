# pyunigen: bindings to the UniGen almost uniform sampler

This directory provides Python bindings to UniGen on the C++ level,
i.e. when importing pycryptosat, the CryptoMiniSat solver becomes part of the
Python process itself.

## Installing

```
pip install pyunigen
```

## Compiling
If you don't want to use the pip package, you can compile it as:

```
apt-get install python-dev
cd python
git clone https://github.com/msoos/cryptominisat
git clone https://github.com/meelgroup/arjun
git clone https://github.com/meelgroup/approxmc
cd ..
python -m build
```
You will then find the files under "dist/".

## Usage

The `pyunigen` module has one object, `Sampler` that has two functions
`sample` and `add_clause`.

The funcion `add_clause()` takes an iterable list of literals such as
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

You can give the following arguments to `Counter`:
* `seed` -- sets the random seed
* `verbosity` -- sets the verbosity of the system (default = 0)
* `epsilon` -- Tolerance parameter, i.e. sets how approximate the returned count is. Default = 0.8
* `delta` -- Confidence parameter, i.e. sets how probabilistically correct the returned count is. Default = 0.20

