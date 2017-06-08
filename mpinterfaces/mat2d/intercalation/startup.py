from __future__ import print_function, division, unicode_literals

import operator

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure

from mpinterfaces.mat2d.intercalation.analysis import get_interstitial_sites


__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


def inject_ions(structure, ion, atomic_fraction):
    """
    Adds ions to a percentage of interstitial sites into a structure
    that results in an at% less than or equal to the specified
    atomic_fraction. Starts by filling interstitial sites with the
    largest voronoi radius, and then works downward.

    Args:
        structure (Structure): Pymatgen Structure object to
            intercalate into.
        ion (str): name of atom to intercalate, e.g. 'Li', or 'Mg'.
        atomic_fraction (int): This fraction of the final intercalated
            structure will be intercalated atoms. Must be < 1.0.

    Returns:
        structure. Includes intercalated atoms.
    """

    specie = Element(ion)

    # If the structure isn't big enough to accomodate such a small
    # atomic fraction, multiply it in the x direction.
    n_ions = 1.
    dimensions = [2, 1, 1]
    while not n_ions / (structure.num_sites+n_ions) <= atomic_fraction:
        structure.make_supercell(dimensions)
        dimensions = [dimensions[1], dimensions[0], 1]

    if structure.num_sites * atomic_fraction > 3:
        print("The algorithm is working, but may take several minutes "
              "due to the relatively large number of ions to "
              "intercalate.")

    interstitial_sites = get_interstitial_sites(structure)["tetrahedral"]

    while n_ions / (structure.num_sites + 1) <= atomic_fraction:
        at_p = int(round(n_ions*100. / (structure.num_sites + 1), 0))
        try:
            structure.append(species=specie, coords=interstitial_sites[0][0],
                             validate_proximity=True,
                             properties={'velocities': [0.0, 0.0, 0.0]},
                             coords_are_cartesian=True)
            interstitial_sites =get_interstitial_sites(structure)["tetrahedral"]
            n_ions += 1

            print("Currently at ~{} at %".format(at_p))

        except ValueError:
            # The largest site is too close to another atom, so your
            # structure is already full!
            raise ValueError("The atomic fraction specified exceeds the "
                             "number of reasonably distant interstitial "
                             "sites in the structure. Please choose a "
                             "smaller atomic fraction and try again.")

    return structure
