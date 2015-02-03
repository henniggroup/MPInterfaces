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
from pymatgen.matproj.rest import MPRester

from mpinterfaces import Interface, get_struct_from_mp


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
    print '\ntransformation matrix : ', t_mat
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
    print 'original u, v', [u1, v1]
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
    print 'reduced u, v', [u, v]
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
