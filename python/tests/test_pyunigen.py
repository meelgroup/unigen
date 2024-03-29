from pycryptosat import Solver
from pyapproxmc import Counter
from pyunigen import Sampler


def test():
    clauses = [ [1, -2, -3], [-1, 3, -2], [2, -1, -3], [1, 3, 2], [4, -5, -6], [-4, -5, 6], 
                [5, -4, -6], [4, 5, 6], [1, -7, -3], [-1, 3, -7], [7, -1, -3], [1, 3, 7], 
                [-8, 4, -6], [-8, -4, 6], [8, -4, -6], [4, 8, 6], [9, 7, 8], [-4, -1, 9, 8], 
                [-4, -1, 9, 7], [-8, 4, -9, -7], [4, -9, -7], [-8, 1, -9, -7], [1, -9, -7], 
                [-8, -9, -7], [-8, 4, -9], [-8, 4, 1, -9], [4, 1, -9, -7], [-8, 1, -9], [10, -1], 
                [10, -3], [1, 3, -10], [10, -11, -12], [-11, -10, 12], [-12, 11, -10], [10, 11, 12], 
                [11, -9, -13], [-11, 9, -13], [-11, -9, 13], [9, 11, 13], [-4, 14], [-6, 14], 
                [4, 6, -14], [-15, 14, -16], [15, -14, -16], [-15, -14, 16], [15, 14, 16], 
                [17, 13, 16], [17, -13, -16], [-17, -13, 16], [-17, 13, -16], [18, -19], [18, -20], 
                [18, -21], [18, -22], [-17, 18], [18, -23], [19, 20, 21, 22, 17, 23, -18], 
                [-24, 25], [-26, -18, 25], [26, 24, -25], [24, 18, -25], [25, -27], [27, -25], 
                [-4, 28], [28, -29], [-15, 28], [28, -30], [28, -31], [-6, 28], [4, 29, 15, 30, 31, 6, -28], 
                [-32, 33], [-34, -28, 33], [32, 34, -33], [32, 28, -33], [33, -35], [35, -33], [36, -37], 
                [-1, 36], [-3, 36], [36, -38], [36, -39], [-12, 36], [37, 1, 3, 38, 39, 12, -36], 
                [-40, 41], [-42, -36, 41], [40, 42, -41], [40, 36, -41], [41, -43], [43, -41], 
                [19, -44], [44, -19], [21, -45], [45, -21], [17, -46], [46, -17], [23, -47], [47, -23], 
                [20, -48], [48, -20], [22, -49], [49, -22], [26, -50], [50, -26], [24, -51], [51, -24], 
                [52, -27], [27, -52], [4, -53], [53, -4], [6, -54], [54, -6], [15, -55], [55, -15], 
                [29, -56], [56, -29], [30, -57], [57, -30], [31, -58], [58, -31], [34, -59], [59, -34], 
                [32, -60], [60, -32], [61, -35], [35, -61], [1, -62], [62, -1], [3, -63], [63, -3], 
                [12, -64], [64, -12], [37, -65], [65, -37], [38, -66], [66, -38], [39, -67], [67, -39], 
                [42, -68], [68, -42], [40, -69], [69, -40], [70, -43], [43, -70], [-11, -16, -71], 
                [9, -16, -71], [-11, 9, -71], [11, 16, 71], [-9, 16, 71], [11, -9, 71], [-4, 72], 
                [-15, 72], [-6, 72], [4, 15, 6, -72], [-29, 72, -73], [29, -72, -73], [-29, -72, 73], 
                [29, 72, 73], [-1, 74], [-3, 74], [-12, 74], [1, 3, 12, -74], [-37, 74, -75], [37, -74, -75], 
                [-37, -74, 75], [37, 74, 75], [73, -75, -76], [71, -75, -76], [71, 73, -76], [-73, 75, 76], 
                [-71, 75, 76], [-71, -73, 76], [-37, 77], [-1, 77], [-3, 77], [-12, 77], [37, 1, 3, 12, -77], 
                [-38, 77, -78], [38, -77, -78], [-38, -77, 78], [38, 77, 78], [-76, 78, -79], [76, -78, -79], 
                [-76, -78, 79], [76, 78, 79], [-4, 80], [-29, 80], [-15, 80], [-6, 80], [4, 29, 15, 6, -80], 
                [-30, 80, -81], [30, -80, -81], [-30, -80, 81], [30, 80, 81], [73, -81, -82], [-73, 81, -82], 
                [82, -73, -81], [73, 81, 82], [20, 79, 82], [20, -79, -82], [-20, -79, 82], [-20, 79, -82], 
                [4, 19, 1, -83], [4, -1, -19, -83], [1, -4, -19, -83], [19, -4, -1, -83], [-4, -1, -19, 83], 
                [19, 1, -4, 83], [4, 19, -1, 83], [4, 1, -19, 83], [-83], [-78, -82, -84], [76, -82, -84], 
                [76, -78, -84], [78, 82, 84], [-76, 82, 84], [-76, 78, 84], [-37, 85], [-1, 85], [-3, 85], 
                [-38, 85], [-12, 85], [37, 1, 3, 38, 12, -85], [-39, 85, -86], [39, -85, -86], [-39, -85, 86], 
                [39, 85, 86], [-4, 87], [-29, 87], [-15, 87], [-30, 87], [-6, 87], [4, 29, 15, 30, 6, -87], 
                [-31, 87, -88], [31, -87, -88], [-31, -87, 88], [31, 87, 88], [-73, 89], [-81, 89], 
                [73, 81, -89], [-88, 89, -90], [88, -89, -90], [-88, -89, 90], [88, 89, 90], [-86, -90, -91],
                [84, -90, -91], [84, -86, -91], [86, 90, 91], [-84, 90, 91], [-84, 86, 91], [-42, 36, -92], 
                [42, -36, -92], [-42, -36, 92], [42, 36, 92], [-34, 28, -93], [34, -28, -93], [-34, -28, 93], 
                [34, 28, 93], [-73, 94], [-81, 94], [-88, 94], [73, 81, 88, -94], [-93, 94, -95], [93, -94, -95], 
                [-93, -94, 95], [93, 94, 95], [-92, 95, -96], [91, 95, -96], [91, -92, -96], [92, -95, 96],
                [-91, -95, 96], [-91, 92, 96], [-37, 97], [-1, 97], [-3, 97], [-38, 97], [-39, 97], [-12, 97],
                [-42, 97], [37, 1, 3, 38, 39, 12, 42, -97], [40, 96, 97, -98], [-40, -96, 97, -98], [40, -96, -97, -98],
                [-40, 96, -97, -98], [-40, -96, -97, 98], [40, 96, -97, 98], [-40, 96, 97, 98], [40, -96, 97, 98],
                [93, 99], [94, 99], [73, 81, 88, 93, 99], [-73, -93, -94, -99], [-81, -93, -94, -99], 
                [-93, -94, -99], [-88, -93, -94, -99], [-4, 100], [-34, 100], [-29, 100], [-15, 100], 
                [-30, 100], [-31, 100], [-6, 100], [4, 34, 29, 15, 30, 31, 6, -100], [32, 99, 100, -101], 
                [-32, -99, 100, -101], [32, -99, -100, -101], [-32, 99, -100, -101], [-32, -99, -100, 101], 
                [32, 99, -100, 101], [-32, 99, 100, 101], [32, -99, 100, 101], [24, 98, 101], [24, -98, -101], 
                [-24, -98, 101], [-24, 98, -101], [-71, 75, -102], [71, -75, -102], [-71, -75, 102], [71, 75, 102], 
                [23, 73, 102, -103], [-23, 73, -102, -103], [-23, -73, 102, -103], [23, -73, -102, -103], 
                [-23, -73, -102, 103], [23, -73, 102, 103], [23, 73, -102, 103], [-23, 73, 102, 103], 
                [-103], [-91, 92, -104], [91, -92, -104], [-91, -92, 104], [91, 92, 104], [26, 95, 104, -105], 
                [-26, 95, -104, -105], [-26, -95, 104, -105], [26, -95, -104, -105], [-26, -95, -104, 105], 
                [26, -95, 104, 105], [26, 95, -104, 105], [-26, 95, 104, 105], [-105], [4, 106], [1, 106], 
                [-4, -1, -106], [-7, 106, -107], [7, -106, -107], [-7, -106, 107], [7, 106, 107], [21, 8, 107, -108], 
                [-8, -21, 107, -108], [8, -21, -107, -108], [-8, 21, -107, -108], [-8, -21, -107, 108], 
                [21, 8, -107, 108], [-8, 21, 107, 108], [8, -21, 107, 108], [-108], [-84, 86, -109], 
                [84, -86, -109], [-84, -86, 109], [84, 86, 109], [22, 90, 109], [22, -90, -109], [-22, 90, -109], [-22, -90, 109]]

    sampling_set = {64, 65, 66, 67, 68, 69, 70, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}

    # Get solution count
    counter = Counter(seed=2157, epsilon=0.5, delta=0.15)
    for clause in clauses:
        counter.add_clause(clause)
    cell_count, hash_count = counter.count()

    # Sample solutions and check that they are valid
    for i in range(1000):
        if i % 100 == 99:
            print("Running test %d" % i)
        # Initialize sampler
        sampler = Sampler(seed=i)
        for clause in clauses:
            sampler.add_clause(clause)
        
        # Generate and extract a sample
        raw_sample = sampler.sample(num=1,sampling_set=sampling_set, cell_hash_count=(cell_count, hash_count))
        sample = raw_sample[2][0]

        # Add sampled values as clauses and check that the formula is still satisfiable
        solver = Solver()
        for clause in clauses:
            solver.add_clause(clause)
        for lit in sample:
            solver.add_clause([lit])
        if solver.solve()[0] == False:
            print("Sample is incorrect: ", sample)
            assert False

if __name__ == '__main__':
    test()
