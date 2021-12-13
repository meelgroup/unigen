from pyapproxmc import Counter
from pyunigen import Sampler

NUM_SAMPLES = 1000
SAMPLE_TOL = 0.20

def minimal_test():
    counter = Counter(seed=2157, epsilon=0.5, delta=0.15)
    counter.add_clause(list(range(1,100)))
    cell_count, hash_count = counter.count()

    sampler = Sampler(seed=2157)
    sampler.add_clause(list(range(1,100)))
    sampler.sample(cell_count, hash_count)

def sampling_set_test():
    counter = Counter(seed=2157, epsilon=0.5, delta=0.15, sampling_set=range(1,10))
    counter.add_clause(list(range(1,100)))
    cell_count, hash_count = counter.count()

    sampler = Sampler(seed=2157, sampling_set=range(1,10))
    sampler.add_clause(list(range(1,100)))
    sampler.sample(cell_count, hash_count)

def real_example_test():
    counter = Counter(seed=120, epsilon=0.01, delta=0.001)
    sampler = Sampler(seed=120, kappa=0.15)

    # Create formula with 10 variables, in which 9 or 10
    # must be true. Should have roughly 50-50 (true-false)
    # split on variables 1-8 and 66-33 split for variables 9, 10.
    counter.add_clause([9,10])
    sampler.add_clause([9,10])

    cell_count, hash_count = counter.count()

    samples = []

    # Sample multiple times this time.
    for i in range(NUM_SAMPLES):
        new_sample = sampler.sample(cell_count, hash_count)
        samples.append([1 if val > 0 else 0 for val in new_sample])

    counts = [sum(sample)/NUM_SAMPLES for sample in zip(*samples)]

    for i in range(8):
        assert 0.5 - SAMPLE_TOL <= counts[i] <= 0.5 + SAMPLE_TOL

    for i in range(8,10):
        assert 0.66 - SAMPLE_TOL <= counts[i] <= 0.66 + SAMPLE_TOL

if __name__ == '__main__':
    minimal_test()
    sampling_set_test()
    real_example_test()
