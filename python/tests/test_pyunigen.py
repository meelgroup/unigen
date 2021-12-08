from pyapproxmc import Counter
from pyunigen import Sampler

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
    counter = Counter(seed=120, epsilon=0.8, delta=0.2, sampling_set=list(range(1,21)))
    sampler = Sampler(seed=120, sampling_set=range(1,21))

    with open("test_1.cnf") as test_cnf:
        # Pop sampling set and metadata lines
        lines = test_cnf.readlines()[2:]

        # Add clauses to counter
        for line in lines:
            literals = [int(i) for i in line.split()[:-1]]
            counter.add_clause(literals)
            sampler.add_clause(literals)

    cell_count, hash_count = counter.count()

    # Sammple multiple times this time.
    for i in range(50):
        sampler.sample(cell_count, hash_count)

if __name__ == '__main__':
    minimal_test()
    sampling_set_test()
    real_example_test()
