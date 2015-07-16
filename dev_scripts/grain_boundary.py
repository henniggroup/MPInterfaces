from mpinterfaces.interface import Interface
import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mpinterfaces import get_struct_from_mp
#initial bulk structure to start with


def get_grain_boundary_interface(structure=None, hkl_pair= {'hkl': [[1,0,0],[1,1,0]],\
                                 'thickness':[10,10]}, twist = 0, tilt = 0, separation=0):
    """
    Args:
        structure: pymatgen structure to create grain boundary in
        hkl_pair:  dict of {'hkl':thickness}
        twist:     twist in degrees
        tilt:      tilt in degrees
    """


    structure = get_struct_from_mp(structure, MAPI_KEY="")
    sa = SpacegroupAnalyzer(structure)
    structure_conventional = sa.get_conventional_standard_structure()
    structure = structure_conventional.copy()
    structure.sort()

    #creation of lower part of grain boundary
    lower= Interface(structure,\
    hkl = hkl_pair['hkl'][0],
    min_thick = hkl_pair['thickness'][0],
    min_vac = separation+hkl_pair['thickness'][1],
    primitive = False, from_ase = True, center_slab=False)

    lower.to(fmt="poscar", filename="POSCAR_lower.vasp")

    #creation of upper part of grain boundary
    upper= Interface(structure,\
    hkl = hkl_pair['hkl'][1],
    min_thick = hkl_pair['thickness'][1],
    min_vac = 0,
    primitive = False, from_ase = True)

    #find top atoms reference of lower part of gb
    substrate_top_z = np.max(np.array([site.coords for site in lower])[:,2])

    # define twist and tilt vectors
    twist_shift_normal = lower.lattice.matrix[2,:]/\
                         np.linalg.norm(lower.lattice.matrix[2,:])
    tilt_normal = upper.lattice.matrix[1,:]/\
                  np.linalg.norm(upper.lattice.matrix[2,:])

    #define twist operation SymmOp object
    twist_op = SymmOp.from_axis_angle_and_translation(axis= twist_shift_normal,\
                angle=twist, angle_in_radians=False,translation_vec=(0, 0, 0))
    #define tilt operation SymmOp object
    tilt_op = SymmOp.from_axis_angle_and_translation(axis= tilt_normal,\
                angle=tilt, angle_in_radians=False,translation_vec=(0, 0, 0))
    upper.apply_operation(twist_op)
    upper.to(fmt="poscar", filename="POSCAR_upper.vasp")
    upper.apply_operation(tilt_op)

    #define shift separation along twist vector normal to upper plane
    shift = -1*twist_shift_normal/np.linalg.norm(twist_shift_normal) * separation
    #define origin to shift w.r.t top of the lower grain 
    origin = np.array([0,0, substrate_top_z])
    #shift sites in upper 
    for site in upper:
        new_coords = site.coords - origin  +  shift
        lower.append(site.specie, new_coords, coords_are_cartesian=True)
    return lower 

if __name__=='__main__':
    #testing different separations for creating the grain boundary
    for i in [11,12,13,14,15,16,17,18,19,10]:
        strt= get_grain_boundary_interface(structure='CdTe', \
               hkl_pair= {'hkl': [[1,0,0],[1,1,0]],'thickness':[10,10]},\
               twist = 0, tilt = 10, separation=i)
    
        strt.to(fmt='poscar', filename='test_gb/POSCAR_final'+str(i)+'.vasp')
