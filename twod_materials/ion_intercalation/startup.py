from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.defects.point_defects import (
    Interstitial, ValenceIonicRadiusEvaluator
    )


def inject_ions(ion, atomic_fraction):
    """
    Adds ions to a percentage of interstitial sites into the POSCAR
    that results in an at% less than or equal to the specified
    atomic_fraction.

    args:
          specie (str): name of ion to intercalate
          atomic_fraction (int): < 1.0
    """

    specie = Element(ion)
    structure = Structure.from_file('POSCAR')

    n_ions = 1.
    while not n_ions / structure.num_sites <= atomic_fraction:
        structure.make_supercell([2, 1, 1])

    evaluator = ValenceIonicRadiusEvaluator(structure)
    interstitial = Interstitial(structure, radii=evaluator.radii,
                                valences=evaluator.valences)

    interstitial_sites = [site._fcoords for site in interstitial._defect_sites]

    while n_ions / (structure.num_sites + 1) <= atomic_fraction:
        try:
            structure.append(species=specie,
                             coords=interstitial_sites[int(n_ions) - 1],
                             validate_proximity=True)
        except IndexError:
            raise ValueError('The atomic_fraction specified exceeds the '
                             'number of available interstitial sites in this '
                             'structure. Please choose a smaller '
                             'atomic_fraction.')
        n_ions += 1

    return structure
