## Import des modules n√©cessaires ##

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

## Donnees du probl√®me ##

M_eau = 0.01801528 #kg/mol
M_air = 0.028965 #kg/mol
v_in = 20 #m/s
T_0 = 283 #K
P_0 = 10**5 #Pa
R = 8.314 #KJ/mol
g = 9.81 #m/s^2
B = (v_in**2)/2+(7/2*(R*T_0)/M_air)

def courant_init(v_max, dim, exp = False):
    mat_courant = np.zeros((dim, dim))
    v_j = v_max / dim
    for j in range(dim):
        mat_courant[dim - 1 - j,0] = v_j
        mat_courant[dim - 1 - j,dim - 1] = v_j
        v_j += v_max / dim
    for i in range(dim):
        mat_courant[0,i] = v_max
    return mat_courant

## M√©thodes math√©matiques & num√©riques ##

def laplacien_relax(matrice, montagne):
    dim = matrice.shape[0]
    x = sp.symbols('x')
    for b in range(25):
        for i in range(1, dim - 1):
            a=dim-1-i
            for j in range(1, dim - 1):
                    if a<=montagne.subs(x, j):
                        if i<= dim//2:
                            matrice[a,j] = 0
                    else:
                        matrice[i,j] = 1/4 * (matrice[i - 1,j] + matrice[i + 1,j] + matrice[i,j - 1] + matrice[i,j + 1]-4*matrice[i,j])
    matrice[j+1]=0
    return matrice

def v_x(mat) :
    a = mat.shape[0]
    mat_v_x = np.zeros((a-2,a-2))
    for i in range(0,mat_v_x.shape[0]) :
        for j in range(0,mat_v_x.shape[0]) :
            mat_v_x[i,j] = (mat[i,j]-mat[i,j+2])/2
    return mat_v_x

def v_z(mat) :
    a = mat.shape[0]
    mat_v_z = np.zeros((a - 2,a - 2))
    for i in range(0, mat_v_z.shape[0]):
        for j in range(0, mat_v_z.shape[0]):
            mat_v_z[j,i] = -(mat[j,i] - mat[j,i + 2]) / 2
    return mat_v_z

## Condition thermodynamique de la formation du nuage ##

def enthalpie(v_1,v_2,z) :
    return B-((v_1**2+v_2**2)/2+g*z)

def Temp(h) :
    return ((h*M_air)/R)*2/7

def Pr(T) :
    return T**(7/2)*P_0*T_0**(-7/2)

def P_sat(T) :
    return ((((T-273)/40)+1)**2)*10**3

def P_e(r,P) :
    return (r*P)/((M_eau/M_air)+r)

def mat_P_sat(v_1,v_2) :
    mat = np.zeros((v_1.shape[0],v_1.shape[0]))
    for i in range(v_1.shape[0]) :
        for j in range(v_2.shape[0]):
            mat[i,j] = P_sat(Temp(enthalpie(v_1[i,j],v_2[i,j],v_1.shape[0]-j-1)))
    return mat
    
def mat_P_e(r,v_1,v_2) :
    mat = np.zeros((v_1.shape[0],v_1.shape[0]))
    for i in range(v_1.shape[0]) :
        for j in range(v_1.shape[0]) :
            mat[i,j] = P_e(r,Pr(Temp(enthalpie(v_1[i,j],v_2[i,j],v_1.shape[0]-j-1))))
    return mat

def pititnuage(r,v_1,v_2) :
    mat_partielle = mat_P_e(r,v_1,v_2)
    mat_saturante = mat_P_sat(v_1,v_2)
    return mat_partielle > mat_saturante

##Programme de simulation ##

# def initialisation :
#     return None

# def sim :
#     pass

## Fonctions d'affichage ##

def ec_mont(ecran, montagne):
    for j in range(ecran.shape[1]):
        for i in range(ecran.shape[0]):
            if i < montagne.subs('x', j):
                ecran[ecran.shape[0] - 1 - i,j] = 1
    return ecran

def ecran(dim,montagne,cd_nuage): #Fonction qui créée une matrice pour l'affichage
    ecran = np.zeros((dim-2,dim-2))
    ec_mont(ecran, montagne)
    for  i in range(dim-2):
        a=dim-1-i
        for j in range(dim-2):
            if cd_nuage[i][j] == True :
                ecran[i][j] = 2
            if cd_nuage[i][j] == False :
                ecran[i][j] = 1
            if a < montagne.subs("x",j) :
                    ecran[i][j] = 0
    return ecran

def affiche_matrice(ecran):
    lignes, colonnes = ecran.shape
    cases_nulles=ecran==0
    couleurs=np.zeros((lignes, colonnes, 3))
    couleurs[cases_nulles]=[0,0,0]
    couleurs[~cases_nulles]=[0.5,0.7,1.0]        
    ecran=plt.imshow(couleurs)
    ecran=plt.axis('off')
    return ecran




## Param√®tres modifiables par l'utilisateur ##

dim = 100 #Dimension de la matrice (pr√©cision de la discr√©tisation)
v = 10/3.6 #R√©partition des vitesses selon l'altitude (Fonction Sympy en z)
montagne = sp.sympify(f"-(x/4-{dim}/8)**2+2*{dim}/3") #Forme de la montagne selon la position horizontale (Fonction Sympy en x)
r = 1/100 #Taux d'humidit√© de l'air en fonction de l'altitude et du temps (Fonction Sympy en z et t)

## Programme principal ##

mat_courant = courant_init(v,dim)
mat_courant = laplacien_relax(mat_courant, montagne)
v_1 = v_x(mat_courant)
v_2 = v_z(mat_courant)
mP_e = mat_P_e(r, v_1, v_2)
mP_sat = mat_P_sat(v_1, v_2)
cd = pititnuage(r,v_1,v_2)
mat_ecran = ecran(dim, montagne, cd)
print("fin")
