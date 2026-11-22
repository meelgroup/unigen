import itertools

import pytest
from pyunigen import Sampler

# var 7 is in no clause: Arjun takes it out of the CNF into the count's multiplier
CLAUSES = [[1, 2], [-1, 3], [-3, 4, 5], [-2, -5], [4, -6], [7, -7]]
NVARS = 7


def solutions(proj):
    sols = set()
    for bits in itertools.product([False, True], repeat=NVARS):
        if all(any(bits[abs(lit) - 1] == (lit > 0) for lit in cl) for cl in CLAUSES):
            sols.add(tuple(v if bits[v - 1] else -v for v in proj))
    return sols


def make_sampler(arjun):
    s = Sampler(seed=3, arjun=arjun)
    for cl in CLAUSES:
        s.add_clause(cl)
    return s


@pytest.mark.parametrize("arjun", [True, False])
@pytest.mark.parametrize("proj", [[1, 2, 3, 4], [1, 4, 6, 7]])
def test_samples_are_solutions(arjun, proj):
    cells, hashes, samples = make_sampler(arjun).sample(num=50, sampling_set=proj)
    sols = solutions(proj)
    assert len(samples) == 50
    assert all(tuple(s) in sols for s in samples)
    assert cells * 2**hashes == len(sols)


@pytest.mark.parametrize("arjun", [True, False])
def test_no_sampling_set(arjun):
    cells, hashes, samples = make_sampler(arjun).sample(num=20)
    sols = solutions(range(1, NVARS + 1))
    assert len(samples) == 20
    assert all(tuple(s) in sols for s in samples)
    assert cells * 2**hashes == len(sols)


@pytest.mark.parametrize("arjun", [True, False])
def test_unsat(arjun):
    s = Sampler(arjun=arjun)
    s.add_clause([1])
    s.add_clause([-1])
    assert s.sample(num=5) == (0, 0, [])
