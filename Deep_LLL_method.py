def Gram_Schmidt(n, B):
    GS = Matrix(RR, n) # Gram-Schmidt vector
    U = Matrix(RR, n) # Gram-Schmidt coefficients

    for i in range(n):
        GS[i] = B[i]
        U[i, i] = 1.0
        for j in range(i):
            U[i, j] = B[i].inner_product(GS[j]) / (GS[j].norm() ** 2)
            GS[i] = GS[i] - (U[i, j] * GS[j])
    return GS, U


def DeepLLL(n, B, delta):
    BB = vector(RR, n)
    GS, U = Gram_Schmidt(n, B)
    for i in range(n):
        BB[i] = GS[i].norm() ** 2
    #step2
    k = 1
    while k < n:
        for j in range(k)[::-1]:
            #step 5
            if abs(U[k, j]) > 0.50:
                q = round(U[k, j])
                B[k] -= q * B[j]
                for l in range(j + 1):
                    U[k, l] -= q * U[j, l]
        C = B[k].norm() ** 2
        i = 0
        while i < k:
            if C >= delta * BB[i]:
                C = C - U[k, i] ** 2 * BB[i]
                i = i + 1
            else:
                v = B[k]
                for j in range(i+1, k + 1)[::-1]:
                    B[j] = B[j-1]
                B[i] = v
                GS, U = Gram_Schmidt(n, B)
                for i in range(n):
                    BB[i] = GS[i].norm() ** 2
                k = max(i, 1) - 1
        k = k+1
    return B

n = 6
B = Matrix(ZZ, 6)
B[0, 0] = 9172
B[1, 0] = 4226
B[2, 0] = 5559
B[3, 0] = 4292
B[4, 0] = 3923
B[5, 0] = 4946
for i in range(1, n):
    B[i, i] = 1
print('B = ', B)
delta = float(input('delta'))
print(DeepLLL(n, B, delta))
print(' ')

n1 = 8
B1 = Matrix(ZZ, 8)
B1[0, 0] = 7871
B1[1, 0] = 6339
B1[2, 0] = 4436
B1[3, 0] = 3265
B1[4, 0] = 223
B1[5, 0] = 6803
B1[6, 0] = 1417
B1[7, 0] = 7787
for i in range(1, n1):
    B1[i, i] = 1
print('B1 = ', B1)
print(DeepLLL(n1, B1, delta))

'''
output

B =  [9172    0    0    0    0    0]
[4226    1    0    0    0    0]
[5559    0    1    0    0    0]
[4292    0    0    1    0    0]
[3923    0    0    0    1    0]
[4946    0    0    0    0    1]
delta0.75
[ 0  1  0  0  0  1]
[-3  0 -2  1  3 -1]
[ 7 -1  1  0  2  0]
[ 2 -3 -1  3  1  4]
[-4 -2  2 -2  6  2]
[-1  0 -6 -2 -1  0]
 
B1 =  [7871    0    0    0    0    0    0    0]
[6339    1    0    0    0    0    0    0]
[4436    0    1    0    0    0    0    0]
[3265    0    0    1    0    0    0    0]
[ 223    0    0    0    1    0    0    0]
[6803    0    0    0    0    1    0    0]
[1417    0    0    0    0    0    1    0]
[7787    0    0    0    0    0    0    1]
[-2  0  1  1  0  0  0 -2]
[ 0  0 -1  2  1  1 -1 -2]
[ 1  1  1 -1  2  0  0  1]
[ 0 -1  1  1  1  0 -1  2]
[ 2  3  0  1  0  0  1  1]
[-3  2  1  0 -1  1  0  1]
[ 0 -1  2  1  1  1  2 -1]
[-2  0 -2 -1  2 -1  2  1]

'''