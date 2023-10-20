import numpy as np
import sympy as sp

def courant_init(v_max, dim):
    mat_courant = np.zeros((dim, dim))
    v_j = v_max / dim
    for j in range(dim):
        mat_courant[dim - 1 - j][0] = int(v_j)
        mat_courant[dim - 1 - j][dim - 1] = int(v_j)
        v_j += v_max / dim
    for i in range(dim):
        mat_courant[0][i] = v_max
    return mat_courant

def laplacien_relax(matrice, montagne):
    dim = matrice.shape[0]
    x = sp.symbols('x')
    for i in range(1, dim - 1):
        for j in range(1, dim - 1):
            if i > montagne.subs(x, j):
                matrice[dim - 1 - i][j] = round(1/4 * (matrice[i - 1][j] + matrice[i + 1][j] + matrice[i][j - 1] + matrice[i][j + 1]), 2)
            else:
                matrice[i][j] = 0
    return matrice

def ec_mont(ecran, montagne):
    for j in range(ecran.shape[1]):
        for i in range(ecran.shape[0]):
            if i < montagne.subs('x', j):
                ecran[ecran.shape[0] - 1 - i][j] = 1
    return ecran

v_max = 10
dim = 11
mat_courant = courant_init(v_max, dim)
montagne_expr = sp.sympify("-(x-5)**2+5")
ecran = np.zeros((dim, dim))

ecran = ec_mont(ecran, montagne_expr)
result = laplacien_relax(mat_courant, montagne_expr)

print(result)
