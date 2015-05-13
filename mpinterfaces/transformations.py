from __future__ import division, unicode_literals, print_function

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
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def get_trans_matrices(n):
    """
    yields a list of 2x2 transformation matrices for the given supercell
    size n
    """
    def factors(n0):
        for i in range(1, n0+1):
            if n0 % i == 0:
                yield i
    for i in factors(n):
        m = n // i
        yield [ [[i, j], [0, m]] for j in range(m) ]

def get_uv(ab, t_mat):
    """
    return u and v, the supercell lattice vectors obtained through the
    transformation matrix
    """
    print('\ntransformation matrix : ', t_mat)
    u = np.array(ab[0])*t_mat[0][0] + np.array(ab[1])*t_mat[0][1] 
    v = np.array(ab[1])*t_mat[1][1] 
    return [u, v]

def get_reduced_uv(uv):
    """
    returns reduced lattice vectors
    """
    is_not_reduced =  True
    u = np.array(uv[0])
    v = np.array(uv[1])
    u1 = copy(u)
    v1 = copy(v)
    print('original u, v', [u1, v1])
    while is_not_reduced:
        if np.dot(u, v) <0:
            v = -v
        #print 'v after dot', v
        if np.linalg.norm(u) > np.linalg.norm(v):
            #print 'norm u and v', np.linalg.norm(u), np.linalg.norm(v)
            u1 = copy(v)
            v1 = copy(u)
        elif np.linalg.norm(v) > np.linalg.norm(u+v):
            #print 'norm v and v+u', np.linalg.norm(v), np.linalg.norm(v+u)
            v1 =  v+u
        elif np.linalg.norm(v) > np.linalg.norm(u-v):
            #print 'norm v and u-v', np.linalg.norm(v), np.linalg.norm(u-v)
            v1 = v-u
        else:
            is_not_reduced = False
        u = copy(u1)
        v = copy(v1)
    print('reduced u, v', [u, v])
    return [u, v]

def reduced_supercell_vectors(ab, n):
    """
    returns all possible reduced in-plane lattice vectors for the
    given starting unit cell lattice vectors(ab) and the
    supercell size n
    """
    uv_list = []
    for r_tm in get_trans_matrices(n):
        for tm in r_tm:
            uv = get_uv(ab, tm)
            uv = get_reduced_uv(uv)
            uv_list.append(uv)
    return uv_list

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
    print('rmax1, rmax2', rmax1, rmax2)
    print('a2/a1', area2/area1)
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

def get_matching_lattices(iface1, iface2, max_area = 100,
                          max_mismatch = 0.01, max_angle_diff = 1,
                          r1r2_tol= 0.02):
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
    
    print('area1, area2', area1, area2)
    r_list = get_r_list(area1, area2, max_area, tol=r1r2_tol)
    if not r_list:
        print('r_list is empty. Try increasing the max surface area or/and the other tolerance paramaters')
        sys.exit()
    for r1r2 in r_list:
        #print('r1, r2', r1r2, '\n')
        uv1_list = reduced_supercell_vectors(ab1, r1r2[0])
        uv2_list = reduced_supercell_vectors(ab2, r1r2[1])
        for uv1 in uv1_list:
            for uv2 in uv2_list:                
                u_mismatch = get_mismatch(uv1[0], uv2[0])
                v_mismatch = get_mismatch(uv1[1], uv2[1])
                angle1 = get_angle(uv1[0], uv1[1])
                angle2 = get_angle(uv2[0], uv2[1])
                print('\nu1, u2', np.linalg.norm(uv1[0]), np.linalg.norm(uv2[0]))
                print('v1, v2', np.linalg.norm(uv1[1]), np.linalg.norm(uv2[1]))
                print('angle1, angle2', angle1, angle2)
                print('u and v mismatches', u_mismatch, v_mismatch)
                if  abs(u_mismatch) < max_mismatch and abs(v_mismatch) < max_mismatch:
                    if  abs(angle1 - angle2) < max_angle_diff:
                        print('\n FOUND ONE')
                        print('u mismatch in percentage = ', u_mismatch)
                        print('v mismatch in percentage = ', v_mismatch)
                        print('angle1, angle diff', angle1, abs(angle1 - angle2))
                        print('uv1, uv2', uv1, uv2)
                        return uv1, uv2
