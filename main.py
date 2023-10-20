from DAHU.maths import Expression, Matrice

def courant_init(v_max,dim) : #Applique les conditions initiales dans la matrice de courant
    mat_courant = Matrice(dim, dim)
    v_j = v_max/mat_courant.dim()[0]
    for j in range(mat_courant.dim()[0]) :
        mat_courant[mat_courant.dim()[1]-1-j][0] = int(v_j)
        mat_courant[mat_courant.dim()[1]-1-j][mat_courant.dim()[1]-1] = int(v_j)
        v_j += v_max/mat_courant.dim()[0]
    for i in range(mat_courant.dim()[1]) :
        mat_courant[0][i] = v_max
    return mat_courant

def laplacien_relax(matrice,montagne) : #Utilisation de l'Algorithme de Jacobi et prenant en compte la forme de la montagne
    for i in range(1,matrice.nbl-1) :
        for j in range(1,matrice.nbc-1) :
            if i > montagne.eval("x",j) :
                matrice[matrice.nbl-1-i][j] = round(1/4*(matrice[i-1][j]+matrice[i+1][j]+matrice[i][j-1]+matrice[i][j+1]),2)
            else :
                matrice[i][j] = 0
    return matrice


mat_courant = courant_init(10,11)
montagne = Expression("-(x-5)^2+5",["x"])
ecran = Matrice(11,11)

def ec_mont(ecran,montagne) : #Mets la montagne dans la mtrice a partir de son expression
    for j in range(ecran.nbc) :
        for i in range(ecran.nbl) :
            if i < montagne.eval("x",j) :
                ecran[ecran.nbl-1-i][j] = 1
    return ecran



print(laplacien_relax(mat_courant,montagne))
