'''
Lagurange基底簡約アルゴリズム
'''


def Lagurange(B):
    if B[0].norm() > B[1].norm():
        v = B[0]; B[0] = B[1]; B[1] = v

    while(B[0].norm() < B[1].norm()):
        q = -round(B[0].inner_product(B[1])/B[0].norm()**2)
        v = B[1] + q*B[0]
        B[1] = B[0]; B[0] = v

    v = B[0]; B[0] = B[1]; B[1] = v

    return B

print('基底簡約')
B1 = Matrix(ZZ, [[-56, 43], [95, -73]])
B2 = Matrix(ZZ, [[-393, 279, 132], [410, -291, -512]])
B3 = Matrix(ZZ, [[230, -651, 609, -366], [301, -852, 797, -479]])
print('B1 = ', Lagurange(B1))
print('B2 = ', Lagurange(B2))
print('B3 = ', Lagurange(B3))

'''
output
基底簡約
B1 =  [ 1  1]
[-1  2]
B2 =  [  17  -12 -380]
[-393  279  132]
B3 =  [-1 -3 -2 -1]
[ 2 -3  5 -2]
'''