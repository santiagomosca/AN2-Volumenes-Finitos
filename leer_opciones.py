#! /usr/bin/env python3
"""
Pequeño script para pasar opciones
al programa 'tp2_main.py'
"""

import argparse
import sys

# --- #

def opciones():
    descripcion="TP 2 Análisis Numérico II.\n\
    'tp2_main.py': programa para resolver la ecuación\n\
    de Laplace en un perfil IPN 600."
    
    parser = argparse.ArgumentParser(description=descripcion)
    
    parser.add_argument('-t', '--tamElem',
            help='Longitud característica de los puntos de la malla.\
                    afecta el tamaño de los elementos. Por defecto\
                    es <5>.',
            required=False, action='store', type=float, default=5)
    
    parser.add_argument('-m', '--funcMallado',
            help='Función de mallado. Opción <1> genera mallas\
            no estructuradas. El solver resuelve la ecuación de\
            Laplace correctamente. La opción <2> genera mallas\
            estructuradas o no estructuradas; las últimas actualmente\
            no son resueltas adecuadamente.',
            required=True, action='store', type=int)
    
    parser.add_argument('-e', '--estruct',
            help='Si la función de mallado es opción <2>, elegir si la malla\
            es estructurada (S) o no (opción por defecto, N).',
            required=False, action='store', type=str, default='N')
    
    
    args = parser.parse_args()
    
    tm = args.tamElem
    func_malla = args.funcMallado
    if func_malla not in (1,2):
        print('Función de mallado mal especificada. Debe ser <1> o <2>.')
        sys.exit(1)
    tipo_malla = args.estruct

    return tm, func_malla, tipo_malla