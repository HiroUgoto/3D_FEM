import numpy as np

def tetra_gauss(n):    
    if n == 1:
        # Degree of precision: 1
        xi = np.array([
            [1/4, 1/4, 1/4, 1/4]
        ])
        w = np.array([1/6])

    elif n == 4:
        # Degree of precision: 2
        sqrt5 = np.sqrt(5.0)
        a = (5.0 + 3.0 * sqrt5) / 20.0
        b = (5.0 - sqrt5) / 20.0
        
        # (a, b, b, b) の巡回
        xi = np.array([
            [a, b, b, b],
            [b, a, b, b],
            [b, b, a, b],
            [b, b, b, a]
        ])
        w = np.full(4, 1/24)

    elif n == 5:
        # Degree of precision: 3
        xi = np.array([
            [1/4, 1/4, 1/4, 1/4],
            [1/2, 1/6, 1/6, 1/6],
            [1/6, 1/2, 1/6, 1/6],
            [1/6, 1/6, 1/2, 1/6],
            [1/6, 1/6, 1/6, 1/2]
        ])
        w = np.array([-4/30, 9/120, 9/120, 9/120, 9/120])

    elif n == 11:
        # Degree of precision: 4
        sqrt5_14 = np.sqrt(5.0 / 14.0)
        a = (1.0 + sqrt5_14) / 4.0
        b = (1.0 - sqrt5_14) / 4.0
                
        xi = np.zeros((11, 4))
        w = np.zeros(11)
        
        # Index 0
        xi[0] = [1/4, 1/4, 1/4, 1/4]
        w[0] = -74/5625
        
        # Index 1-4 (Type 2)
        val_11 = 11/14
        val_1 = 1/14
        xi[1:5] = np.array([
            [val_11, val_1, val_1, val_1],
            [val_1, val_11, val_1, val_1],
            [val_1, val_1, val_11, val_1],
            [val_1, val_1, val_1, val_11]
        ])
        w[1:5] = 343/45000
        
        # Index 5-10 (Type 3: permutations of a,a,b,b)
        xi[5:11] = np.array([
            [a, a, b, b],
            [a, b, a, b],
            [a, b, b, a],
            [b, a, a, b],
            [b, a, b, a],
            [b, b, a, a]
        ])
        w[5:11] = 56/2250
        
    else:
        raise ValueError(f"Unsupported number of points n={n}. Available: 1, 4, 5, 11")

    return xi[:,0:3], w