#%%
import numpy as np
import sys
from leer_GMSH import xnod_from_msh as extraer_nodos
from leer_GMSH import elem_from_msh as extraer_elementos
from obtener_grupos_fisicos import grupos_fisicos as obtener_grupos_fisicos
from obtener_grupos_fisicos import obtener_nodos as obtener_nodos_grupo
#%%
def leer_malla(filename):
    """
    Función que extrae, del archivo .msh, las coordenadas xyz 
    de los nodos globales de la malla ,y los nodos globales asociados 
    a cada elemento triangular m.
    """
    xyz_nod = extraer_nodos(filename+'.msh')
    ele_tagnod = extraer_elementos(filename+'.msh')    
    datos = {
        "nodos_xyz": xyz_nod,
        "elementos": ele_tagnod,
    }
    return datos
#%
def obtener_nodos_elemento(tag_elem, datos):
    ''' 
    Función que extrae las coordenadas de los nodos locales (1), (2) y (3) 
    del elemento triangular m
    '''
    nlocal = datos["elementos"][tag_elem]
    n1 = datos["nodos_xyz"][nlocal[0]-1]
    n2 = datos["nodos_xyz"][nlocal[1]-1]
    n3 = datos["nodos_xyz"][nlocal[2]-1]
    n = np.array([n1, n2, n3])
    return n, nlocal
#%
def matriz_jacobiano(tag_elem, tipo_elem_triangular, datos):
    '''
    Función que calcula el Jacobiano a partir de un elemento normalizado triangular.
    '''
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
    '''
    Función que calcula las integrales d{Psi}_Eps y d{Psi}_Etta asociados al nodo local k.
    (Calcula la longitud del contorno del sub-dominio del nodo local k)    
    Fila #k corresponde al nodo local #k con k entre 1 y 3
    '''
    if tipo_elemento_triangular == 1:
        # # Triangulo 1 (del mallador GMSH)
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
    '''
    Función que obtiene los elementos de la matriz Jacobiana y el Jacobiano 
    asociados al elemento normalizado
    '''
    matriz, jacobiano = matriz_jacobiano(tag_elem, tipo_elem_triangular, datos)
    L11, L12 = matriz[0,:]
    L21, L22 = matriz[1,:]
        
    alfa1 = +(L22*L21+L11*L12)/jacobiano
    alfa2 = +(L22**2+L12*L12)/jacobiano
    alfa3 = -(L11**2+L21*L21)/jacobiano
    alfa4 = -(L12*L11+L22*L21)/jacobiano
    return alfa1, alfa2, alfa3, alfa4
#%
def obtener_contribuciones_elemento(tag_elem, tipo_elem_triangular,  datos):
    '''
    Función que obtiene las contribuciones del elemento m a los nodos locales k
    (Contribución del nodo local k' al nodo k , ambos del elemento m)
    en una matriz de rigidez local A [k x k]
    
    '''
    # Se obtienen los datos del elemento normalizado elegido
    psi_eps_etta, dN_deps, dN_detta = elem_normalizado(tipo_elem_triangular)

    # Se obtienen los coeficientes alfas
    alfa1, alfa2, alfa3, alfa4 = obtener_coef_elemento(tag_elem, tipo_elem_triangular, datos)

    k= 3 #nodos locales del elemento triangular
    A = np.zeros((k,k))
    for i in range(k):
        # Contribuciones del elemento m al nodo local k
        alfa_eps = (alfa1*psi_eps_etta[i,0])*dN_deps + (alfa2*psi_eps_etta[i,1])*dN_deps
        alfa_etta = (alfa3*psi_eps_etta[i,0])*dN_detta + (alfa4*psi_eps_etta[i,1])*dN_detta
        A[i,:] = alfa_eps+alfa_etta # Matriz que guarda las contribuciones al nodo local k en sus filas
    return A
#%
def obtener_matriz_global(M, N, nodosk, A):
    '''
    Función que obtiene la Matriz de rigidez global M [N x N], a partir los nodos globales
    asociados a un elemento m
    'nodosk': es el diccionario que asocia los nodos locales de un elemento m a los nodos globales de la malla.    
    '''
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
    '''
    Función que aplica las condiciones de contorno de Dirichlet a la Matriz global.
    También devuelve el termino independiente o fuente del sistema.
    'dict_fisico': contiene el index y los nombres de los grupos fisicos.
    'dict_nodos': contiene el mismo index del grupo fisico y los nodos globales asociados.
    'dict_condiciones': contiene los nombres de los grupos físicos y las condiciones a aplicar.
    
    OBS: Los nombres de los grupos fisicos en la malla y los dados en el diccionario condiciones, 
    deben coincidir.
    '''
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
#%
def escribir_resultado(array_resultados, filename):
    """Función que escribe al archivo .msh los resultados de la simulación.
    Anexa al final del archivo el escalar temperatura. Se puede abrir
    directamente con el programa 'Gmsh'. El formato de .msh es 4.1.
    Información sobre el formato para escribir en
    https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

    La cantidad de nodos se obtiene del tamaño del array de resultados.
    Debido a que la indexación de Gmsh comienza con '1', debe añadirse
    la unidad a la numeración de numpy.
    Los resultados de temperatura corresponde a los valores del array.
    """

    with open(filename+'.msh', 'a') as archivo:
        archivo.write('$NodeData\n')
        archivo.write('1\n')    # 1 string tag
        archivo.write('"Temperatura"\n')  # the name of the view
        archivo.write('1\n')    # 1 real tag
        archivo.write('0.0\n')  # the time value
        archivo.write('3\n')    # 3 integer tags:
        archivo.write('0\n')    # the time step
        archivo.write('1\n')    # 1-component (scalar) field
        archivo.write(f"{array_resultados.size}\n") # N associated nodal values
        for nodo, valor in np.ndenumerate(array_resultados):
            archivo.write("{n} {v}\n".format(n=nodo[0]+1, v=valor))
        archivo.write('$EndNodeData\n')
