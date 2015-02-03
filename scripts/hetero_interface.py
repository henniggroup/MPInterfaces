"""
Compute the reduced matching lattice vectors for heterostructure
interfaces as described in the paper by Zur and McGill:
Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
"""

import sys
from math import sqrt
from copy import copy

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from mpinterfaces import Interface, get_struct_from_mp
from mpinterfaces.transformations import reduced_supercell_vectors


def get_r_list(area1, area2, max_area, tol=0.02):
    """
    returns a list of r1 and r2 values that satisfies:
    r1/r2 = area2/area1 with the constraints:
    r1 <= Area_max/area1 and r2 <= Area_max/area2
    r1 and r2 corresponds to the supercell sizes of the 2 interfaces
    that align them     
    """
    r_list = []
    rmax1 = int(max_area/area1)
    rmax2 = int(max_area/area2)
    print 'rmax1, rmax2', rmax1, rmax2
    print 'a2/a1', area2/area1
    for r1 in range(1, rmax1+1):
        for r2 in range(1, rmax2+1):
            if abs(float(r1)/float(r2) - area2/area1) < tol:
                r_list.append([r1, r2])
    return r_list

def get_mismatch(a, b):
    """
    percentage mistmatch between the lattice vectors a and b
    """
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(b)/np.linalg.norm(a) - 1

def get_angle(a, b):
    """
    angle between lattice vectors a and b in degrees
    """
    a = np.array(a)
    b = np.array(b)
    return np.arccos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)) * 180 / np.pi

def get_matching_lattices(iface1, iface2, max_area=200,
                          max_mismatch=0.01, max_angle_diff = 1):
    """
    computes a list of matching reduced lattice vectors that satify
    the max_area, max_mismatch and max_anglele_diff criteria
    """
    if iface1 is None and iface2 is None:
        #test : the numbers from the paper
        a1 = 5.653
        a2 = 6.481
        #for 100 plane
        ab1 = [ [0, a1/2, -a1/2], [0, a1/2,a1/2]]
        ab2 = [ [0, a2/2, -a2/2], [0, a2/2,a2/2]]
        area1 = a1**2 / 2
        area2 = a2**2 / 2
    
        #for 110 plane
        ab1 = [ [a1/2,-a1/2,0], [0, 0,a1]]
        ab2 = [ [a2/2,-a2/2,0], [0, 0,a2]]    
        area1 = a1**2 / sqrt(2)
        area2 = a2**2 / sqrt(2)
    
        #for 111 surface
        #ab1 = [ [a1/2, 0, a1/2], [a1/2, a1/2, 0]]
        #ab2 = [ [a2/2, 0, a2/2], [a2/2, a2/2, 0]]
        #area1 = a1**2 * sqrt(3)/4 #/ 2 /sqrt(2)
        #area2 = a2**2 * sqrt(3)/4 #/ 2 / sqrt(2)
    else:
        area1 = iface1.surface_area
        area2 = iface2.surface_area
    
        #a, b vectors that define the surface    
        ab1 = [ iface1.lattice.matrix[0,:] , iface1.lattice.matrix[1,:] ]
        ab2 = [ iface2.lattice.matrix[0,:] , iface2.lattice.matrix[1,:] ]
    
    print 'area1, area2', area1, area2
    r_list = get_r_list(area1, area2, max_area, tol=0.02)
    print r_list
    #sys.exit()
    for r1r2 in r_list:
        print 'r1, r2', r1r2, '\n'
        uv1_list = reduced_supercell_vectors(ab1, r1r2[0])
        uv2_list = reduced_supercell_vectors(ab2, r1r2[1])
        for uv1 in uv1_list:
            for uv2 in uv2_list:                
                u_mismatch = get_mismatch(uv1[0], uv2[0])
                v_mismatch = get_mismatch(uv1[1], uv2[1])
                angle1 = get_angle(uv1[0], uv1[1])
                angle2 = get_angle(uv2[0], uv2[1])
                print '\nu1, u2', np.linalg.norm(uv1[0]), np.linalg.norm(uv2[0])
                print 'v1, v2', np.linalg.norm(uv1[1]), np.linalg.norm(uv2[1])
                print 'angle1, angle2', angle1, angle2
                print 'u and v mismatches', u_mismatch, v_mismatch
                if  abs(u_mismatch) < max_mismatch and abs(v_mismatch) < max_mismatch:
                    if  abs(angle1 - angle2) < max_angle_diff:
                        print '\n FOUND ONE'
                        print 'u mismatch in percentage = ', u_mismatch
                        print 'v mismatch in percentage = ', v_mismatch
                        print 'angle1, angle diff', angle1, abs(angle1 - angle2)
                        print 'uv1, uv2', uv1, uv2
                        return uv1, uv2

        
if __name__ == '__main__':
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    
    #structure from materials project, use your own key
    #gaas = get_struct_from_mp('GaAs', MAPI_KEY="dwvz2XCFUEI9fJiR")
    #sa_gaas = SpacegroupAnalyzer(gaas)
    #gaas_cvn = sa_gaas.get_conventional_standard_structure()
    #gaas_cvn.to(fmt='poscar', filename='POSCAR_GaAs.vasp')
    
    #cdte = get_struct_from_mp('CdTe', MAPI_KEY="dwvz2XCFUEI9fJiR")
    #sa_cdte = SpacegroupAnalyzer(cdte)
    #cdte_cvn = sa_cdte.get_conventional_standard_structure()
    #cdte_cvn.to(fmt='poscar', filename='POSCAR_CdTe.vasp')    

    #initial conventional bulk unit cells
    gaas_cvn = Structure.from_file('POSCAR.mp-2534_GaAs')
    cdte_cvn = Structure.from_file('POSCAR.mp-406_CdTe')    

    #slabs
    iface_gaas = Interface(gaas_cvn, hkl=[1,1,0], min_thick=5, min_vac=25)
    iface_cdte = Interface(cdte_cvn, hkl=[1,1,0], min_thick=5, min_vac=5)

    #matching lattice
    uv_gaas, uv_cdte = get_matching_lattices(iface_gaas, iface_cdte,
                                             max_area=100,
                                             max_mismatch=0.01,
                                             max_angle_diff=1)

    #create supercells
    gaas_tmp = np.array( [ uv_gaas[0][:], uv_gaas[1][:],
                            iface_gaas.lattice.matrix[2,:] ] )        
    cdte_tmp = np.array( [ uv_cdte[0][:], uv_cdte[1][:],
                            iface_cdte.lattice.matrix[2,:] ] )
    gaas_latt = Lattice(gaas_tmp)    
    cdte_latt = Lattice(cdte_tmp)
    _, __, scell = iface_gaas.lattice.find_mapping(gaas_latt,
                                               ltol = 0.01, atol=1)
    iface_gaas.make_supercell(scell)
    _, __, scell = iface_cdte.lattice.find_mapping(cdte_latt,
                                               ltol = 0.01, atol=1)
    iface_cdte.make_supercell(scell)
    
    ###################
    #change the lattice 
    ##################
    print 'lattice mapping'
    lmap_matrix = np.array( [ iface_cdte.lattice.matrix[0,:],
                              iface_cdte.lattice.matrix[1,:],
                            iface_gaas.lattice.matrix[2,:] ] )
    lmap = Lattice(lmap_matrix)    
    #lattice, _, scell = iface_gaas.lattice.find_mapping(lmap,
    #                                           ltol = 0.01, atol=1)
    #print iface_gaas.lattice
    print 'new lattice',lmap
    iface_gaas.to(fmt='poscar', filename='POSCAR_GaAs_iface_1.vasp')
    iface_gaas.modify_lattice(lmap)
    #iface_gaas.make_supercell(scell)
    print 'final lattices'
    print iface_gaas.lattice
    print iface_cdte.lattice    
    iface_gaas.to(fmt='poscar', filename='POSCAR_GaAs_iface_2.vasp')
    iface_cdte.to(fmt='poscar', filename='POSCAR_CdTe_iface2.vasp')

    #merge
    for site in iface_cdte:
        vec = iface_gaas.lattice.matrix[2,:]
        print vec
        print site.coords
        new_coords = site.coords - vec/np.linalg.norm(vec)
        print new_coords
        iface_gaas.append(site.specie, new_coords, coords_are_cartesian=True)

    iface_gaas.to(fmt='poscar', filename='POSCAR_combo.vasp')        
    
