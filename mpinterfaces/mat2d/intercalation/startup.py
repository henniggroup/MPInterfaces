from __future__ import print_function, division, unicode_literals

import operator

from monty.dev import requires

from pymatgen.analysis.defects.point_defects import Interstitial, \
    ValenceIonicRadiusEvaluator
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

try:
    import zeo
    zeo_found = True
except ImportError:
    zeo_found = False


@requires(zeo_found, 'get_voronoi_nodes requires Zeo++ cython extension to be '
          'installed. Please contact developers of Zeo++ to obtain it.')
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
    while not n_ions / structure.num_sites <= atomic_fraction:
        structure.make_supercell([2, 1, 1])

    evaluator = ValenceIonicRadiusEvaluator(structure)
    interstitial = Interstitial(structure, radii=evaluator.radii,
                                valences=evaluator.valences)

    interstitial_sites = [
        (site._fcoords, site.properties.get('voronoi_radius', None))
        for site in interstitial._defect_sites]

    # Sort the interstitial sites by their voronoi radii.
    interstitial_sites.sort(key=operator.itemgetter(1))
    interstitial_sites.reverse()

    i = 0
    while n_ions / (structure.num_sites + 1) <= atomic_fraction:
        try:
            structure.append(species=specie, coords=interstitial_sites[i][0],
                             validate_proximity=True)
            n_ions += 1
            i += 1

            evaluator = ValenceIonicRadiusEvaluator(structure)
            interstitial = Interstitial(structure, radii=evaluator.radii,
                                        valences=evaluator.valences)

            interstitial_sites = [
                (site._fcoords, site.properties.get('voronoi_radius', None))
                for site in interstitial._defect_sites]

            # Sort the interstitial sites by their voronoi radii.
            interstitial_sites.sort(key=operator.itemgetter(1))
            interstitial_sites.reverse()

        except ValueError:
            i += 1

        except IndexError:
            raise ValueError('The atomic_fraction specified exceeds the '
                             'number of available interstitial sites in this '
                             'structure. Please choose a smaller '
                             'atomic_fraction.')

    return structure
