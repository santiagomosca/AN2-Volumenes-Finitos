#%% Imports
import gmsh
import numpy as np
def cuadrado(altura, base, tamaño_malla, filename):  
    #%%
    # filename = "cuadrado"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(filename)

    #%% Perfil cuadrado
    tm = tamaño_malla
    tmr = 0.5

    h = altura
    b = base

    geom = gmsh.model.geo

    P1 = geom.addPoint(0,0,0,tm, tag = 1)
    p1 = (0, P1)
    l_ext = geom.extrude([p1], b, 0, 0)
    p2 = l_ext[0]
    l1 = l_ext[1]
    l_ext = geom.extrude([p2], 0, h, 0)
    p3 = l_ext[0]
    l2 = l_ext[1]
    l_ext = geom.extrude([p3], -b, 0, 0)
    p4 = l_ext[0]
    l3 = l_ext[1]
    l4 = geom.addLine(p4[1], p1[1])

    geom.synchronize()

    contorno = geom.addCurveLoop([l1[1], l2[1], l3[1], l4])
    superficie = geom.addPlaneSurface([contorno])

    geom.synchronize()

    s = gmsh.model.addPhysicalGroup(2, [superficie])
    gmsh.model.setPhysicalName(2, s, "Superficie")

    T1 = gmsh.model.addPhysicalGroup(1, [l4])
    gmsh.model.setPhysicalName(1, T1, "Lateral izquierdo")

    T3 = gmsh.model.addPhysicalGroup(1, [l2[1]])
    gmsh.model.setPhysicalName(1, T3, "Lateral derecho")

    T2 = gmsh.model.addPhysicalGroup(1, [l1[1]])
    gmsh.model.setPhysicalName(1, T2, "Borde inferior")
    
    G1 = gmsh.model.addPhysicalGroup(1, [l3[1]])
    gmsh.model.setPhysicalName(1, G1, "Borde superior")

    geom.synchronize()

    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)

    gmsh.write(filename+".msh")
    gmsh.finalize()
    
    return filename
#%%
def perfil_ipn(largo_base, ancho_base, altura_total, espesor_alma, tamaño_malla, filename):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(filename)

    #%% Creación del perfil
    tm = tamaño_malla
    tmr = 0.5
    h = altura_total
    b = largo_base
    s = espesor_alma
    t = ancho_base
    # h = 600.0
    # b = 215.0
    # s = 21.0
    # t = 33.0

    geom = gmsh.model.geo

    p1 = geom.addPoint(0, 0, 0, tm, tag = 1)
    l_ext1 = geom.extrude([(0, p1)], b, 0, 0)
    p2 = l_ext1[0]
    l1 = l_ext1[1]
    l_ext2 = geom.extrude([p2], 0, t, 0)
    p3 = l_ext2[0]
    l2 = l_ext2[1]
    #l_ext3 = geom.extrude([p3], -(b-2*s)/2, 0, 0)
    l_ext3 = geom.extrude([p3], -(b-s)/2, 0, 0)
    p4 = l_ext3[0]
    l3 = l_ext3[1]
    l_ext4 = geom.extrude([p4], 0, h-2*t, 0)
    p5 = l_ext4[0]
    l4 = l_ext4[1]
    #l_ext5 = geom.extrude([p5], (b-2*s)/2, 0, 0)
    l_ext5 = geom.extrude([p5], (b-s)/2, 0, 0)
    p6 = l_ext5[0]
    l5 = l_ext5[1]
    l_ext6 = geom.extrude([p6], 0, t, 0)
    p7 = l_ext6[0]
    l6 = l_ext6[1]
    l_ext7 = geom.extrude([p7], -b, 0, 0)
    p8 = l_ext7[0]
    l7 = l_ext7[1]
    l_ext8 = geom.extrude([p8], 0, -t, 0)
    p9 = l_ext8[0]
    l8 = l_ext8[1]
    #l_ext9 = geom.extrude([p9], (b-2*s)/2, 0, 0)
    l_ext9 = geom.extrude([p9], (b-s)/2, 0, 0)
    p10 = l_ext9[0]
    l9 = l_ext9[1]
    l_ext10 = geom.extrude([p10], 0, -(h-2*t), 0)
    p11 = l_ext10[0]
    l10 = l_ext10[1]
    #l_ext11 = geom.extrude([p11], -(b-2*s)/2, 0, 0)
    l_ext11 = geom.extrude([p11], -(b-s)/2, 0, 0)
    p12 = l_ext11[0]
    l11 = l_ext11[1]
    l12 = geom.addLine(p12[1], p1)

    geom.synchronize()

    contorno = geom.addCurveLoop([l1[1], l2[1], l3[1], l4[1], l5[1], l6[1], l7[1], l8[1], l9[1], l10[1], l11[1], l12])
    superficie = geom.addPlaneSurface([contorno])

    geom.synchronize()

    s = gmsh.model.addPhysicalGroup(2, [superficie]) 
    gmsh.model.setPhysicalName(2, s, "Superficie") 

    T2 = gmsh.model.addPhysicalGroup(1, [l1[1]])
    gmsh.model.setPhysicalName(1, T2, "Borde inferior")

    T1 = gmsh.model.addPhysicalGroup(1, [l2[1], l3[1], l4[1], l5[1], l6[1]])
    gmsh.model.setPhysicalName(1, T1, "Lateral derecho")

    G1 = gmsh.model.addPhysicalGroup(1, [l7[1]])
    gmsh.model.setPhysicalName(1, G1, "Borde superior")

    T3 = gmsh.model.addPhysicalGroup(1, [l8[1], l9[1], l10[1], l11[1], l12])
    gmsh.model.setPhysicalName(1, T3, "Lateral izquierdo")

    geom.synchronize()

    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)  # Ver las "caras" de los elementos finitos 2D
    gmsh.option.setNumber('Mesh.Points', 1)        # Ver los nodos de la malla

    gmsh.write(filename+".msh")
    gmsh.finalize()
    return filename

#%%
def estruct_perfil_ipn(largo_base, ancho_base, altura_total, espesor_alma, tamaño_malla, filename, estructurada='N', refinado=0):
    """
    Creación de malla para perfil IPN.
    La malla puede ser estructurada con elementos triangulares.
    En tal caso la opción es 'S'. Por defecto es no ('N').
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(filename)

    #--- Parámetros del perfil ---#

    h = altura_total
    b = largo_base
    s = espesor_alma
    t = ancho_base


    #--- Parámetros de la malla ---#
 
    tm = tamaño_malla

    # Cantidad de nodos en curvas regulares
    if estructurada=='S':
        el_extr_ala=[t/tm]
        el_extr_alma=[(h-t*2)/tm]
    else:
        el_extr_ala=[]
        el_extr_alma=[]


    #--- Creación de geometría y mallado ---#

    # Ala inferior: puntos y curvas
    p1 = gmsh.model.geo.addPoint(0, 0, 0, tm)
    p2 = gmsh.model.geo.addPoint((b/2-s/2), 0, 0, tm)
    p3 = gmsh.model.geo.addPoint((b/2+s/2), 0, 0, tm)
    p4 = gmsh.model.geo.addPoint(b, 0, 0, tm)
 
    l1 = gmsh.model.geo.addLine(p1,p2)
    l2 = gmsh.model.geo.addLine(p2,p3)
    l3 = gmsh.model.geo.addLine(p3,p4)

    # Ala inferior: extrusión de superficie
    s1 = gmsh.model.geo.extrude([(1,l1)], 0, t, 0,
            numElements=el_extr_ala)  # Extremo izq del ala
    s2 = gmsh.model.geo.extrude([(1,l2)], 0, t, 0,
            numElements=el_extr_ala)  # Centro del ala
    s3 = gmsh.model.geo.extrude([(1,l3)], 0, t, 0,
            numElements=el_extr_ala)  # Extremo der del ala
 
    # Ala superior: puntos y curvas
    p5 = gmsh.model.geo.addPoint(0, (h-t), 0, tm)
    p6 = gmsh.model.geo.addPoint((b/2-s/2), (h-t), 0, tm)
    p7 = gmsh.model.geo.addPoint((b/2+s/2), (h-t), 0, tm)
    p8 = gmsh.model.geo.addPoint(b, (h-t), 0, tm)
 
    l4 = gmsh.model.geo.addLine(p5,p6)
    l5 = gmsh.model.geo.addLine(p6,p7)
    l6 = gmsh.model.geo.addLine(p7,p8)

    # Ala superior: extrusión de superficie
    s4 = gmsh.model.geo.extrude([(1,l4)], 0, t, 0,
            numElements=el_extr_ala)  # Extremo izq del ala
    s5 = gmsh.model.geo.extrude([(1,l5)], 0, t, 0,
            numElements=el_extr_ala)  # Centro del ala
    s6 = gmsh.model.geo.extrude([(1,l6)], 0, t, 0,
            numElements=el_extr_ala)  # Extremo der del ala
 

    # Alma: extrusión
    s7 = gmsh.model.geo.extrude([s2[0]], 0, (h-t*2), 0,
            numElements=el_extr_alma)
 
    # Sincronización de la geometría
    gmsh.model.geo.synchronize()
    
    # Mallado 2D de la geometría
    gmsh.model.mesh.generate(2)
    
    #--- Refinación de la malla ---#
    # Refina la malla dividiendo en 2 los elementos generados previamente.
    # La refinación es igual a 0 por defecto
    for ref in range(refinado):
        gmsh.model.mesh.refine()

    #--- Identificación de líneas y superficies físicas ---#
    T1 = gmsh.model.addPhysicalGroup(1, [s3[2][1], s3[0][1], s7[2][1], l6, s6[2][1]])
    gmsh.model.setPhysicalName(1, T1, "Lateral derecho")

    T2 = gmsh.model.addPhysicalGroup(1, [l1,l2,l3])
    gmsh.model.setPhysicalName(1, T2, "Borde inferior")
 
    T3 = gmsh.model.addPhysicalGroup(1, [s1[3][1],s1[0][1], s7[3][1], l4, s4[3][1]])
    gmsh.model.setPhysicalName(1, T3, "Lateral izquierdo")
    
    G1 = gmsh.model.addPhysicalGroup(1, [s4[0][1],s5[0][1],s6[0][1]])
    gmsh.model.setPhysicalName(1, G1, "Borde superior")
    
    sups = []
    for sup in gmsh.model.getEntities(2):
        sups.append(sup[1])
    
    SUPS= gmsh.model.addPhysicalGroup(2, sups)
    gmsh.model.setPhysicalName(2, SUPS, "Superficie") 

    # --- Opciones de visualización --- #

    gmsh.option.setNumber('Mesh.SurfaceFaces', 0)  # Ver las "caras" de los elementos finitos 2D
    gmsh.option.setNumber('Mesh.Points', 0)        # Ver los nodos de la malla


    # --- Escritura a archivo y cierre --- #

    gmsh.write(filename+".msh")
    gmsh.finalize()
    
    # --- Fin de función --- #
   
    return filename
