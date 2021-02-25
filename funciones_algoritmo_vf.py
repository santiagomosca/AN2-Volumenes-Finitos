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
def obtener_nodos_elemento(tag_elem, datos):
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
    coord_nodo_local, nodosk = obtener_nodos_elemento(tag_elem, datos)
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
    dN_deps = np.array([0, -1, 1])
    dN_detta = np.array([-1, 1, 0])
    # # Triangulo 1 
    # psi_eps_etta = np.array([[-.5, .5], [.0, -.5], [.5, .0]])
    return psi_eps_etta, dN_deps, dN_detta
#%
def obtener_coef_elemento(tag_elem, datos):
    matriz, jacobiano = matriz_jacobiano(tag_elem, datos)
    lx23, lx12 = matriz[0,:]
    ly23, ly12 = matriz[1,:]
    
    alfa1 = +(ly12*ly23+lx23*ly23)/jacobiano
    alfa2 = +(ly12**2+ly23*lx12)/jacobiano
    alfa3 = -(lx23**2+lx12*ly23)/jacobiano
    alfa4 = -(lx12*lx23+ly12*lx12)/jacobiano
    return alfa1, alfa2, alfa3, alfa4
#%
def obtener_contribuciones_elemento(tag_elem, datos):
    psi_eps_etta, dN_deps, dN_detta = elem_normalizado()

    alfa1, alfa2, alfa3, alfa4 = obtener_coef_elemento(tag_elem, datos)
    
    k= 3 #nodos locales
    A = np.zeros((k,k))
    for i in range(k):
        # Contribuciones del elemento m al nodo local i
        alfa_eps = (alfa1*psi_eps_etta[i,0])*dN_deps + (alfa2*psi_eps_etta[i,1])*dN_detta
        alfa_etta = (alfa3*psi_eps_etta[i,0])*dN_deps + (alfa4*psi_eps_etta[i,1])*dN_detta
        A = alfa_eps+alfa_etta
    return A
#%
def obtener_matriz_global(M, N, nodosk, A):
    nodosk = [(nodosk[k]-1) for k in range(len(nodosk))]
    k = 3 #nodos locales
    for it in range(k):
        n = nodosk[it]
        # a_mn = A[it]
        a_mn = np.zeros(N)
        a_mn[nodosk] = A[it]
        M[:,n-1] += a_mn
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
    I = np.eye(num_nodos)
    matriz[U_cond.nonzero(),:] = I[U_cond.nonzero()]
    fuente = U_cond
    return matriz, fuente