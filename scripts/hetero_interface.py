"""
Compute the reduced matching lattice vectors for heterostructure interfaces
as described in the paper by Zur and McGill:
Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
"""

import numpy as np
from mpinterfaces.interface import Interface
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.matproj.rest import MPRester
from copy import copy
import sys

MAPI_KEY="dwvz2XCFUEI9fJiR"

def get_struct_from_mp(formula):
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        print "\nnumber of structures matching the chemical formula "+formula+" = ", len(data)
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            return structure

def get_r_list(area1, area2, max_area, tol=0.01):
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
    for r1 in range(1, rmax1+1):
        for r2 in range(1, rmax2+1):
            if float(r1)/float(r2) - area2/area1 < tol:
                r_list.append([r1, r2])
    return r_list

def get_im(n):
    """
    returns list of i and m such that i*m = n
    where n is the supercell size
    """
    m_list = []
    i_list = []
    for i in range(1, n+1):
        for m in range(1, n+1):
                if i*m == n:
                    m_list.append(m)
                    i_list.append(i)
    return (i_list, m_list)

def get_trans_matrices(n):
    """
    returns a list of 2x2 transformation matrices for the given supercell
    size n
    """
    i_list, m_list = get_im(n)
    trans_matrices = []
    for k in range(len(i_list)):
        for j in range(m_list[k]):
            t_mat = [ [i_list[k], j], [0, m_list[k]] ]
            trans_matrices.append(t_mat)
    return trans_matrices
    

def get_uv(ab, t_mat):
    """
    return u and v, the supercell lattice vectors obtained through the 
    transformation matrix
    """
    u = np.array(ab[0])*t_mat[0][0] + np.array(ab[1])*t_mat[0][1] 
    v = np.array(ab[1])*t_mat[1][1] 
    return [u, v]

def get_reduced_uv(uv):
    """
    reduced lattice vectors
    """
    is_not_reduced =  True
    u = np.array(uv[0])
    v = np.array(uv[1])
    u1 = copy(u)
    v1 = copy(v)
    #print 'original u, v', [u1, v1]
    while is_not_reduced:
        if np.dot(u, v) <0:
            v = -v
        if np.linalg.norm(u) > np.linalg.norm(v):
            u1 = copy(v)
            v1 = copy(u)
        elif np.linalg.norm(v) > np.linalg.norm(u+v):
            u1 = copy(u)
            v1 =  v+u
        elif np.linalg.norm(v) > np.linalg.norm(u-v):
            u1 = copy(u)
            v1 = v-u
        else:
            is_not_reduced = False
        u = copy(u1)
        v = copy(v1)
    #print 'reduced u, v', [u, v]
    return [u, v]

def get_mismatch(a, b):
    """
    mistmatch between the lattice vectors a and b in percentage
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

def get_matching_lattices(iface1, iface2, max_area=50, max_mismatch=1, max_angle_diff = 1):
    """
    computes a list of matching reduced lattice vectors that satify the max_area
    max_mismatch and max_anglele_diff criteria
    """
    area1 = iface1.surface_area
    ab1 = [ iface1.lattice.matrix[0,:] , iface1.lattice.matrix[1,:] ] 
    area2 = iface2.surface_area
    ab2 = [ iface2.lattice.matrix[0,:] , iface2.lattice.matrix[1,:] ] 
    print 'area1, area2', area1, area2
    r_list = get_r_list(area1, area2, max_area)
    for r1r2 in r_list:
        print 'r1, r2', r1r2, '\n'
        r1_t_mats = get_trans_matrices(r1r2[0])
        r2_t_mats = get_trans_matrices(r1r2[1])
        uv1_list = []
        uv2_list = []
        for r1_tm in r1_t_mats:
            uv1 = get_uv(ab1, r1_tm)
            uv1 = get_reduced_uv(uv1)
            uv1_list.append(uv1)
        for r2_tm in r2_t_mats:
            uv2 = get_uv(ab2, r2_tm)
            uv2 = get_reduced_uv(uv2)
            uv2_list.append(uv2)
        for uv1 in uv1_list:
            for uv2 in uv2_list:
                u_mismatch = get_mismatch(uv1[0], uv2[0])
                v_mismatch = get_mismatch(uv1[1], uv2[1])
                if  abs(u_mismatch) < max_mismatch and abs(v_mismatch) < max_mismatch:
                    angle1 = get_angle(uv1[0], uv1[1])
                    angle2 = get_angle(uv2[0], uv2[1])
                    if  abs(angle1 - angle2) < max_angle_diff:
                        print 'u mismatch in percentage = ', u_mismatch
                        print 'v mismatch in percentage = ', v_mismatch
                        print 'angle1, angle diff', angle1, abs(angle1 - angle2)
                        print 'uv1, uv2', uv1, uv2
                        sys.exit()

        
if __name__ == '__main__':
    strt1 = get_struct_from_mp('GaAs')
    strt2 = get_struct_from_mp('CdTe')

    iface1 = Interface(strt1, hkl=[1,1,1], min_thick=10, min_vac=5)
    iface2 = Interface(strt2, hkl=[1,1,1], min_thick=10, min_vac=5)

    get_matching_lattices(iface1, iface2, max_area=100, max_mismatch=0.1, max_angle_diff=1)

#    print get_r_list(21, 15.98, 200, tol=0.01)
#    print get_trans_matrices(4)
