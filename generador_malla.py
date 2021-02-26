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
    l_ext3 = geom.extrude([p3], -(b-2*s)/2, 0, 0)
    p4 = l_ext3[0]
    l3 = l_ext3[1]
    l_ext4 = geom.extrude([p4], 0, h-2*t, 0)
    p5 = l_ext4[0]
    l4 = l_ext4[1]
    l_ext5 = geom.extrude([p5], (b-2*s)/2, 0, 0)
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
    l_ext9 = geom.extrude([p9], (b-2*s)/2, 0, 0)
    p10 = l_ext9[0]
    l9 = l_ext9[1]
    l_ext10 = geom.extrude([p10], 0, -(h-2*t), 0)
    p11 = l_ext10[0]
    l10 = l_ext10[1]
    l_ext11 = geom.extrude([p11], -(b-2*s)/2, 0, 0)
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