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
def matriz_jacobiano(tag_elem, tipo_elem_triangular, datos):
    # Calcula el Jacobiano a partir de un elemento normalizado triangular
    coord_nodo_local, nodosk = obtener_nodos_elemento(tag_elem, datos)
    x = [coord_nodo_local[k,0] for k in range(len(nodosk))]
    y = [coord_nodo_local[k,1] for k in range(len(nodosk))]
    psi, dN_deps, dN_detta = elem_normalizado(tipo_elem_triangular)
    L11, L12 = np.dot(dN_deps, x), np.dot(dN_detta, x)
    L21, L22 = np.dot(dN_deps, y), np.dot(dN_detta, y)
    matriz = np.array([[L11, L12], [L21, L22]])
    jacobiano = L11*L22 - L12*L21
    return matriz, jacobiano
#%
def elem_normalizado(tipo_elemento_triangular):
    # Calcula las integrales d{Psi}_Eps y d{Psi}_Etta del nodo
    # Fila #k corresponde al nodo local #k con k entre 1 y 3
    if tipo_elemento_triangular == 1:
        # # Triangulo 1 
        psi_eps_etta = np.array([[-.5, .5], [.0, -.5], [.5, .0]])
        dN_deps = np.array([-1, 1, 0])
        dN_detta = np.array([-1, 0, 1])   
    elif tipo_elemento_triangular == 2:
        # Mi triangulo 
        psi_eps_etta = np.array([[.0, .5 ],[-.5, .5],[.5, .0]])
        dN_deps = np.array([0, -1, 1])
        dN_detta = np.array([-1, 1, 0])
    else:
        print("El tipo de elemento no está definido")
        print("")
    return psi_eps_etta, dN_deps, dN_detta
#%
def obtener_coef_elemento(tag_elem, tipo_elem_triangular, datos):
    matriz, jacobiano = matriz_jacobiano(tag_elem, tipo_elem_triangular, datos)
    L11, L12 = matriz[0,:]
    L21, L22 = matriz[1,:]
    
    alfa1 = +(L22*L21+L11*L21)/jacobiano
    alfa2 = +(L22**2+L21*L12)/jacobiano
    alfa3 = -(L11**2+L12*L21)/jacobiano
    alfa4 = -(L12*L11+L22*L12)/jacobiano
    return alfa1, alfa2, alfa3, alfa4
#%
def obtener_contribuciones_elemento(tag_elem, tipo_elem_triangular,  datos):
    psi_eps_etta, dN_deps, dN_detta = elem_normalizado(tipo_elem_triangular)

    alfa1, alfa2, alfa3, alfa4 = obtener_coef_elemento(tag_elem, tipo_elem_triangular, datos)
    # alfa1, alfa2, alfa3, alfa4 = obtener_coef_elemento(tag_elem, datos)

    k= 3 #nodos locales
    A = np.zeros((k,k))
    for i in range(k):
        # Contribuciones del elemento m al nodo local i
        alfa_eps = (alfa1*psi_eps_etta[i,0])*dN_deps + (alfa2*psi_eps_etta[i,1])*dN_detta
        alfa_etta = (alfa3*psi_eps_etta[i,0])*dN_deps + (alfa4*psi_eps_etta[i,1])*dN_detta
        A[i,:] = alfa_eps+alfa_etta
    return A
#%
def obtener_matriz_global(M, N, nodosk, A):
    nodosk = [(nodosk[k]-1) for k in range(len(nodosk))]
    k = 3 #nodos locales
    for it in range(k):
        n = nodosk[it]
        # a_mn = A[it]
        a_mn = np.zeros(N)
        a_mn[nodosk] = A[it,:]
        M[n,:] += a_mn
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