def Gram_Schmidt(n, B):
    GS = Matrix(RR, n)
    U = Matrix(RR, n)

    for i in range(n):
        GS[i] = B[i]
        U[i, i] = 1.0
        for j in range(i):
            U[i, j] = B[i].inner_product(GS[j]) / (GS[j].norm() ** 2)
            GS[i] = GS[i] - (U[i, j] * GS[j])
    return GS, U


def check_LLL(n, B, delta):
    # Size-reduce?
    GS, U = Gram_Schmidt(n, B)
    print()
    for i in range(n):
        for j in range(i):
            if abs(U[i, j]) > 0.50:
                return False

    for k in range(1, n):
        if GS[k].norm()**2 < (delta - U[k, k - 1] ** 2) * GS[k - 1].norm() ** 2:
            return False
    return True

    # lovasz condition?


delta = 0.75
n = 4
B = Matrix(ZZ, [[-1, 0, -1, 5],
                [1, -1, -1, -3],
                [-1, -4, 2, 1],
                [3, -1, -3, 2]])

print(check_LLL(n, B, delta))

n = 6
C = Matrix(ZZ, [[0, -17, 8, -7, 9, 6],
                [-9,  9, -7, 14, 13, 1],
                [4, 14, 1, 3, -3, 20],
                [18, -12, 6, 16, 7, 1],
                [19, -3, -20, -11, 11, 13],
                [20, 9, 7, -10, 28, 0]] )
print(check_LLL(n, C, delta))


'''
output

False #B is not reduced 

True  #C is reduced
'''