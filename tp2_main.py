#! /usr/bin/env python3
"""
TP 2 Análisis Numérico II.


Programa para resolver la ecuación de Laplace en un perfil IPN 600.equeño script para pasar opciones
al programa 'tp2_main.py'
"""

import numpy as np
from numpy.__config__ import show
from pandas import DataFrame as show_matrix
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as splinalg
import funciones_algoritmo_vf as funcion
import generador_malla
import leer_opciones
from obtener_grupos_fisicos import grupos_fisicos as obtener_grupos_fisicos

# --- Lectura de opciones --- #

geo_malla, tm, func_malla, tipo_malla, ref_malla = leer_opciones.opciones()

# --- Datos de perfil IPN 600 --- #
altura_total = 600.0
base = 215.0
espesor = 21.0
ancho_base = 33.0

# --- Datos de cuadrado --- #
arg1_cuad = 10
arg2_cuad = 10

# --- Creación de la malla --- #
if geo_malla == 'C':
    filename = generador_malla.cuadrado(arg1_cuad, arg2_cuad, tm, 'cuadrado')

else: # geo_malla=='P':
    filename = 'perfil_ipn'
    
    if func_malla==1:
        # Función no estructurada original
        filename = generador_malla.perfil_ipn(base, ancho_base, altura_total,
                espesor, tm, filename)
    else: # func_malla==2
        # Función para malla estructurada
        filename = generador_malla.estruct_perfil_ipn(base, ancho_base,
                altura_total, espesor, tm, filename, estructurada=tipo_malla, refinado=ref_malla)

datos = funcion.leer_malla(filename)

# --- Condiciones de contorno --- #
condiciones = {
    "Lateral izquierdo": 16,
    "Lateral derecho": 18,
    "Borde inferior": 20,
}

# --- Creación de la matriz --- #
M, N = np.size(datos["elementos"],0), np.size(datos["nodos_xyz"],0)
Matriz_Global = np.zeros((N,N))
print("Cantidad de elementos 'M': {:.0f}".format(M))
print("Cantidad de nodos 'N': {:.0f}".format(N))

# Ciclo for (barrido de elementos)
for m in range(M):
    G_mn, nodos_k = funcion.obtener_nodos_elemento(m, datos)
    A = funcion.obtener_contribuciones_elemento(m,tipo_elem_triangular=1, datos=datos)
    Matriz_Global = funcion.obtener_matriz_global(Matriz_Global, N, nodos_k, A)  

print("La matriz global sin condiciones de contorno:")
show_matrix(Matriz_Global)

# --- Condiciones Iniciales y de contorno --- #
grupo_fisico_dict, nodo_grupo_dict = obtener_grupos_fisicos(filename+'.msh')
Matriz_Global, fuente = funcion.cc_dirichlet(grupo_fisico_dict, nodo_grupo_dict, condiciones, Matriz_Global)

print("La Matriz Global con las condiciones de contorno")
show_matrix(Matriz_Global)
print("")
M = csc_matrix(Matriz_Global)

# --- Resolver el sistema lineal AX = B --- #
# Temperaturas = splinalg.spsolve(M, fuente)
Temperaturas, istop, itn, normr = splinalg.lsqr(M, fuente)[:4]

# --- Escritura al archivo .msh de los resultados, para visualizar con Gmsh --- #
print("Escritura a archivo .msh de las temperaturas")
funcion.escribir_resultado(Temperaturas, filename)
#-----------------------------------------------------------
# %%
