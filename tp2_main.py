#%%
import numpy as np
from numpy.__config__ import show
from pandas import DataFrame as show_matrix
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as splinalg
import funciones_algoritmo_vf as funcion
import generador_malla
from obtener_grupos_fisicos import grupos_fisicos as obtener_grupos_fisicos
# Perfil Ipn
tm = 50

altura_total = 600.0
base = 215.0
espesor = 21.0
ancho_base = 33.0

# Creacion malla
filename = generador_malla.cuadrado(10, 10, 5, 'cuadrado')
# filename = 'perfil_ipn'
# filename = generador_malla.perfil_ipn(base, ancho_base, altura_total, espesor, tm, filename)
datos = funcion.leer_malla(filename)

# Condiciones de contorno
condiciones = {
    "Lateral izquierdo": 16,
    "Borde inferior": 20,
    "Lateral derecho": 18,
}

# Creacion Matriz
M, N = np.size(datos["elementos"],0), np.size(datos["nodos_xyz"],0)
Matriz_Global = np.zeros((N,N))
print("Cantidad de elementos 'M': {:.0f}".format(M))
print("Cantidad de nodos 'N': {:.0f}".format(N))
    # Ciclo for (barrido de elementos)
for m in range(M):
    G_mn, nodos_k = funcion.elem_nodos(m, datos)
    A = funcion.coef_a_mn(m, datos)
    Matriz_Global = funcion.matriz_global(Matriz_Global, N, nodos_k, A)  
print("La matriz global sin condiciones de contorno:")
show_matrix(Matriz_Global)

# Condiciones Iniciales y de contorno  
grupo_fisico_dict, nodo_grupo_dict = obtener_grupos_fisicos(filename+'.msh')
Matriz_Global, fuente = funcion.cc_dirichlet(grupo_fisico_dict, nodo_grupo_dict, condiciones, Matriz_Global)
print("La Matriz Global con las condiciones de contorno")
show_matrix(Matriz_Global)
print("")
M = csc_matrix(Matriz_Global)
#%%
# Resolver el sistema lineal AX = B
Temperaturas = splinalg.spsolve(M, fuente)
# x, istop, itn, normr = lsqr(A, b)[:4]
print("Temperaturas")
print(list(Temperaturas))
#-----------------------------------------------------------
# %%
