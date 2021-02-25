#%%
import numpy as np
import sys
from leer_GMSH import xnod_from_msh as extraer_nodos
from leer_GMSH import elem_from_msh as extraer_elementos
from obtener_grupos_fisicos import grupos_fisicos as obtener_grupos_fisicos
from obtener_grupos_fisicos import obtener_nodos as obtener_nodos_grupo
#%%
def leer_malla(filename):
    xyz_nod = extraer_nodos(filename+'.msh')
    ele_tagnod = extraer_elementos(filename+'.msh')    
    datos = {
        "nodos_xyz": xyz_nod,
        "elementos": ele_tagnod,
    }
    return datos
#%
def elem_nodos(tag_elem, datos):
    # Extrae las coordenadas de los nodos locales (1), (2) y (3) del elemento m
    nlocal = datos["elementos"][tag_elem]
    n1 = datos["nodos_xyz"][nlocal[0]-1]
    n2 = datos["nodos_xyz"][nlocal[1]-1]
    n3 = datos["nodos_xyz"][nlocal[2]-1]
    n = np.array([n1, n2, n3])
    return n, nlocal
#%
def matriz_jacobiano(tag_elem, datos):
    # Calcula el Jacobiano del elemento
    # Mi triangulo
    coord_nodo_local, nodosk = elem_nodos(tag_elem, datos)
    lx23 = coord_nodo_local[2,0]-coord_nodo_local[1,0]
    lx12 = coord_nodo_local[1,0]-coord_nodo_local[0,0]
    ly23 = coord_nodo_local[2,1]-coord_nodo_local[1,1]
    ly12 = coord_nodo_local[1,1]-coord_nodo_local[0,1]
    jacobiano = lx23*ly12 - lx12*ly23
    matriz = np.array([[lx23, lx12],[ly23, ly12]])
    return matriz, jacobiano
#%
def elem_normalizado():
    # Calcula las integrales d{Psi}_Eps y d{Psi}_Etta del nodo
    # Fila #k corresponde al nodo local #k con k entre 1 y 3
    
    # Mi triangulo 
    psi_eps_etta = np.array([[.0, .5 ],[-.5, .5],[.5, .0]])
    # # Triangulo 1 
    # psi_eps_etta = np.array([[-.5, .5], [.0, -.5], [.5, .0]])
    return psi_eps_etta
#%
def coef_alfa(tag_elem, datos):
    matriz, jacobiano = matriz_jacobiano(tag_elem, datos)
    lx23, lx12 = matriz[0,:]
    ly23, ly12 = matriz[1,:]
    
    alfa12 = (lx12**2 + ly12**2)/jacobiano
    alfa23 = (lx23**2 + ly23**2)/jacobiano
    alfa123 = (lx12*lx23 + ly12*ly23)/jacobiano
    return alfa12, alfa23, alfa123
#%
def coef_a_mn(tag_elem, datos):
    psi_eps_etta = elem_normalizado()

    alfa12, alfa23, alfa123 = coef_alfa(tag_elem, datos)
    alfa = np.array([[-alfa23, alfa123], [alfa123, alfa12]])
    
    k= 3 #nodos locales
    A = np.zeros((k,k))
    for i in range(k):
        # Contribuciones del elemento m al nodo local i
        a1_m_i = np.dot(alfa[0,:], psi_eps_etta[i,:])
        a3_m_i = np.dot(alfa[1,:], psi_eps_etta[i,:])
        # Mi triangulo
        if i == 1:
            A[i,0] = a1_m_i
        elif i==2:
            A[i,1] = a1_m_i - a3_m_i
        else:
            A[i,2] = a3_m_i
    return A
#%
def matriz_global(M, N, nodosk, A):
    nodosk = [(nodosk[k]-1) for k in range(len(nodosk))]
    k = 3 #nodos locales
    for it in range(k):
        n = nodosk[it]
        # a_mn = A[it]
        a_mn = np.zeros(N)
        a_mn[nodosk] = A[it,:]
        M[n-1,:] += a_mn
    return M
#%
def cc_dirichlet(dict_fisico, dict_nodos, dict_condiciones, matriz):
    nombres_grupo = [t[1] for t in dict_fisico.values()]
    tag_grupos = list(dict_condiciones.keys())
    num_nodos = np.size([matriz],2)
    U_cond = np.zeros(num_nodos)
    for it in range(len(tag_grupos)):
        if tag_grupos[it] in nombres_grupo:
            nodos_grupo = dict_nodos[list(dict_nodos.keys())[nombres_grupo.index(tag_grupos[it])]]
            nodos_grupo = [int(nodos_grupo[t]-1) for t in range(len(nodos_grupo))]
            u_cond_actual = list(dict_condiciones.values())[it]
            U_cond[nodos_grupo] = u_cond_actual
        else:
            print("\nError en CC_DIRICHLET: No se encuentra el grupo '{}'".format(tag_grupos[it]))
            print("\nLos grupos obtenidos del archivo .msh son:\n")
            print(nombres_grupo)
            print("")
            sys.exit(1)
    fuente = U_cond
    I = np.eye(num_nodos)
    matriz[U_cond.nonzero(),:] = I[U_cond.nonzero()]
    return matriz, fuente