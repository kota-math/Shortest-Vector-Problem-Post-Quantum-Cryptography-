'''
Input: {b1, . . . , bn}: n 次元格子 L の基底，
µi,j (1 ≤ j < i ≤ n), Bi = ∥b∗i∥2(1 ≤ i ≤ n): Gram-Schmidt 情報，R21 ≤ · · · ≤ R2n: 数え上げ上界列
Output: すべての 1 ≤ k ≤ n に対して ∥πn+1−k (v)∥2 ≤ R2k を満たす格子ベクトル v = ∑ni=1 vibi ∈ L の係数ベクトル (v1, . . . , vn) ∈ Zn
1: σ ← (0)(n+1)×n
2: r0 = 0;r1 = 1; · · · ;rn = n
3: ρ1 = · · · = ρn+1 = 0
4: v1 = 1; v2 = · · · = vn = 0
5: c1 = · · · = cn = 0
6: w1 = · · · = wn = 0
7: last nonzero = 1 /∗ vi ̸= 0 となる最大の 1 ≤ i ≤ n ∗/
8: k = 1
1: while true do
2: ρk ← ρk+1 + (vk − ck )2Bk /∗ ρk = ∥πn+1−k (v)∥2 ∗/
3: if ρk ≤ R+1−k then
4: if k = 1 then return (v1, . . . , vn) /∗ v ∈ L の係数ベクトル ∗/
5: k ← k − 1; rk−1 ← max(rk−1,rk )
6: for i = rk downto k + 1 do: σi,k ← σi+1,k + viµi,k ;
7: ck ← −σk+1,k ; vk ← ⌊ck ⌉; wk ← 1
8: else
9: k ← k + 1; if k = n + 1 then return ∅ /∗ v ∈ L が存在しない ∗/
10: rk−1 ← k
11: if k ≥ last nonzero then
12: last nonzero ← k; vk ← vk + 1
13: else
14: if vk > ck then vk ← vk − wk ; else vk ← vk + wk ;
15: wk ← wk + 1
16: end if
17: end if
18: end while

'''

'''
最短ベクトルの数え上げプログラム
基底行列で定まる格子上の最短なベクトルを見つける
'''


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
# ENUM algorithm
def ENUM(n, B, R):
    BB = vector(RR, n)
    GS, U = Gram_Schmidt(n, B)
    for i in range(n):
        BB[i] = GS[i].norm() ** 2
    sigma = Matrix(RR, n+1, n)
    r = vector(ZZ, n+1)
    for i in range(n+1):
        r[i] = 1
    rho = vector(RR, n+1)
    v = vector(ZZ, n)
    v[0] = 1
    c = vector(RR, n)
    w = vector(ZZ, n)
    last_nonzero = 1
    k = 1
    #Main loop
    while(1):
        # Step 2
        rho[k-1] = rho[k] + (v[k-1] - c[k-1]) ** 2 * BB[k-1]
        if rho[k-1] <= R ** 2:
            if k == 1:
                print("Solution founded")
                return v
            # Step 5
            k = k-1
            r[k-1]  = max(r[k-1],r[k])
            # Step 6
            for i in range(k+1, r[k]+1)[::-1]:
                sigma[i-1,k-1] = sigma[i,k-1] + v[i-1]*U[i-1,k-1]
            # Step 7
            c[k-1] = -sigma[k,k-1]
            v[k-1] = round(c[k-1])
            w[k-1] = 1
        else:
            # Step 9
            k = k+1
            if k == n+1 :
                print(0)
                print("No solution.")
                return False
            r[k-1] = k
            if k >= last_nonzero:
                last_nonzero = k
                v[k-1] = v[k-1] + 1
            else:
                # Step 14
                if v[k-1] > c[k-1] :
                    v[k-1] = v[k-1] - w[k-1]
                else:
                    v[k-1] = v[k-1] + w[k-1]
                # Step 15
                w[k-1] = w[k-1] +1
    return v


delta = 0.75
n1 = 4
B1 = Matrix(ZZ, [[18,7,-3,-1],
                [12,-9,-8,15],
                [-9,1,-18,6],
                [-19,-9,0,6],])

n2 = 7
B2 = Matrix(ZZ, [[2,4,18,9,-11,11,2],
                [12,6,19,-26,-7,-22,-20],
                [-13,-16,-23,10,17,24,-20],
                [-16,9,8,8,-9,19,22],
                [19,-8,-23,2,-24,-2,-9],
                [-1,-5,-4,20,-16,3,-11],
                [-27,8,30,26,-3,-13,-5]])
print(B1)
R = 7
v = ENUM(n1, B1 ,R)
if v != False:
    w = vector(ZZ, n1)
    for i in range(n1):
        w += v[i]*B1[i]
    print(w)


print(B2)
R = 16
v = ENUM(n2, B2 ,R)
if v != False:
    w = vector(ZZ, n2)
    for i in range(n2):
        w += v[i]*B2[i]
    print(w)


n3 = 6
B3 = Matrix(ZZ, [[471691,0,0,0,0,0],
                [213169,1,0,0,0,0],
                [376903,0,1,0,0,0],
                [9871,0,0,1,0,0],
                [254689,0,0,0,1,0],
                [150236,0,0,0,0,1]])
B3 = B3.LLL()
print(B3)

R = 0.99 * B3[0].norm()
while(1):
    v = ENUM(n3, B3, R)
    if v != False:
        print(v)
        w = vector(ZZ, n3)
        for i in range(n3):
            w = w + v[i] * B3[i]
        print('w = ', w)
        print(RR(w.norm()))
        R = 0.99 * RR(w.norm())
    else:
        print('End')
        break



n4 = 8
B4 = Matrix(ZZ, [[32121,0,0,0,0,0,0,0],
                [4817,1,0,0,0,0,0,0],
                [2054,0,1,0,0,0,0,0],
                [7023,0,0,1,0,0,0,0],
                [22244,0,0,0,1,0,0,0],
                [8571,0,0,0,0,1,0,0],
                [9357,0,0,0,0,0,1,0],
                [15650,0,0,0,0,0,0,1]])

B4 = B4.LLL()
print(B4)

R = 16
while(1):
    v = ENUM(n4, B4, R)
    if v != False:
        w = vector(ZZ, n4)
        for i in range(n4):
            w += v[i]*B4[i]
        print("w = ", w)
        print("Norm = ", RR(w.norm()))
        R = 0.99*RR(w.norm())
    else:
        break

 '''
 
[ 18   7  -3  -1]
[ 12  -9  -8  15]
[ -9   1 -18   6]
[-19  -9   0   6]
Solution founded
(-1, -2, -3, 5)
[  2   4  18   9 -11  11   2]
[ 12   6  19 -26  -7 -22 -20]
[-13 -16 -23  10  17  24 -20]
[-16   9   8   8  -9  19  22]
[ 19  -8 -23   2 -24  -2  -9]
[ -1  -5  -4  20 -16   3 -11]
[-27   8  30  26  -3 -13  -5]
Solution founded
(3, -6, -3, 4, 13, 2, 0)
 
[ 3  5  3  5 -2  1]
[ 7  1 -4  3  0 -1]
[-3 -3  2  2 -6  3]
[-2 -3 -2  2 -1 -8]
[ 1  7 -3 -7 -3  0]
[ 2 -3  8 -6  0 -6]
Solution founded
(0, 0, 1, 0, 0, 0)
w =  (-3, -3, 2, 2, -6, 3)
8.42614977317636
0
No solution.
End
[ 0 -1  1  0 -1  0  1  1]
[-1  1  1 -1  1  2  1  1]
[ 3 -1 -2 -1  1  0  1 -1]
[-2 -1 -1 -1 -1 -1 -2  2]
[ 2  0 -2  1  1  2 -1  2]
[ 0 -2  0 -1  1  3  0 -2]
[-1 -2 -3 -1 -2  0  2 -1]
[-2 -2  0  3  4  0  3  0]
Solution founded
w =  (0, -1, 1, 0, -1, 0, 1, 1)
Norm =  2.23606797749979
0
No solution.

 
 
 '''