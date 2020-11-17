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

def LLL(n, B, delta):
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

        if BB[k] >= (delta - U[k, k-1]**2) * BB[k-1]:
            k = k + 1
        else:
            v = B[k-1]; B[k-1] = B[k]; B[k] = v
            GS, U = Gram_Schmidt(n, B)
            for i in range(n):
                BB[i] = GS[i].norm() ** 2
            k = max(k-1, 1)

    return B

n = 6
B = Matrix(ZZ, 6)
B[0, 0] = 7612
B[1, 0] = 1162
B[2, 0] = 2835
B[3, 0] = 2860
B[4, 0] = 1129
B[5, 0] = 5712
for i in range(1, n):
    B[i, i] = 1
print(B)
delta = float(input('delta'))
print(LLL(n, B, delta))
print(check_LLL(n, B, delta))

#B1を代入
B1 = Matrix(ZZ, 8)
B1[0, 0] = 32121
B1[1, 0] = 4817
B1[2, 0] = 2054
B1[3, 0] = 7023
B1[4, 0] = 22244
B1[5, 0] = 8571
B1[6, 0] = 9357
B1[7, 0] = 15650
for i in range(1, 8):
    B1[i, i] = 1
print(B1)
print(LLL(8, B1, delta))
print(check_LLL(8, B1, 0.99)) #B1がLLL簡約されているか

'''
output
B =
[7612    0    0    0    0    0]
[1162    1    0    0    0    0]
[2835    0    1    0    0    0]
[2860    0    0    1    0    0]
[1129    0    0    0    1    0]
[5712    0    0    0    0    1]
delta0.99
[ 0  1  1 -3 -1  1]
[ 1 -1 -2 -2  1  2]
[ 1 -2  0  0 -3  1]
[-3  0 -1 -3  0 -2]
[ 5  0 -1 -1  0 -3]
[-2 -3  4 -1  4  1]

True

B1 = 
[32121     0     0     0     0     0     0     0]
[ 4817     1     0     0     0     0     0     0]
[ 2054     0     1     0     0     0     0     0]
[ 7023     0     0     1     0     0     0     0]
[22244     0     0     0     1     0     0     0]
[ 8571     0     0     0     0     1     0     0]
[ 9357     0     0     0     0     0     1     0]
[15650     0     0     0     0     0     0     1]
[ 0 -1  1  0 -1  0  1  1]
[-1  1  1 -1  1  2  1  1]
[ 3 -1 -2 -1  1  0  1 -1]
[-2 -1 -1 -1 -1 -1 -2  2]
[ 2  0 -2  1  1  2 -1  2]
[ 0 -2  0 -1  1  3  0 -2]
[-1 -2 -3 -1 -2  0  2 -1]
[-2 -2  0  3  4  0  3  0]

True
'''