from __future__ import print_function, division, unicode_literals

import os
import operator

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from math import sqrt

from scipy.spatial import Delaunay, ConvexHull

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mpinterfaces.utils import is_converged

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


def sq_dist(p1, p2):
    """
    Calculate the non-square-root distance between two points.

    Args:
        p1, p2: 1x3 point coordinates.
    """
    return (p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2


def pt_btwn(pt1, pt2, r):
    """
    Get the vector of magnitude `r` along the path from pt1 to pt2.
    Args:
        pt1, pt2 (array): points defining the direction of the
            vector to return
        r (float): magnitude of vector to return
    Returns:
        3D vector.
    """
    total_vector = np.subtract(pt2, pt1)
    u = np.array(total_vector/np.linalg.norm(total_vector))
    return np.add(pt1, r*u)


def get_interstitial_sites(structure, octahedra=False, unique=False):
    """
    Use a Delaunay triangulation of all atomic sites in the crystal
    structure to define tetrahedra of open volumes (interstitial
    sites). Each interstitial site is ranked according to the maximum
    radius of an atom that could fit in that site without overlapping
    one of the existing neighboring atoms' radii.

    The default behavior is to stop there, but by setting `octahedra`
    to True, the tetrahedra which share faces are combined to form
    bipyramids (hexahedra) and then points are added to these
    bipyramids to formoctahedra, in order to identify the largest 5-
    and 6-fold coordinated sites as well. This takes a little longer
    since it requires combining tetrahedra.

    Args:
        structure (Structure): Pymatgen Structure object
        octahedra (Boolean): Whether or not to search also for
            octahedral interstitial sites.
        unique (Boolean): Whether or not to enforce that only
            symmetrically inequivalent sites are returned.
            Determining the symmetry-equivalence is usually
            by far the slowest task in the algorithm.
    Returns:
        interstitials (dict): dictionary of the form
            {"tetrahedral": [(coordinates, max_radius), ...],
             "hexahedral": [(coordinates, max_radius), ...],
             "octahedral": [(coordinates, max_radius), ...]}
            storing lists of each interstitial site for both
            coordination types, sorted by largest radius first.
            Coordinates are given as cartesian.
    """

    # Preserve the original structure
    st = structure.copy()

    # Small unit cells make the triangulation unreliable
    n_sites = structure.num_sites
    if n_sites < 4:
        st.make_supercell(3)
    m_0 = st.lattice._matrix

    # Make a 3x3x3 supercell so that the center unit cell
    # is surrounded by its images- i.e. it has no "boundaries",
    # which can erroneously create tetrahedra of infinite volumes.
    st.make_supercell(3)
    m = st.lattice._matrix

    # These are the vertices of only the center cell
    cell_vertices = np.array([
        np.add(np.add(m[0]/3., m[1]/3.), m[2]/3.),
        np.add(np.add(m[0]/1.5, m[1]/3.), m[2]/3.),
        np.add(np.add(m[0]/3., m[1]/1.5), m[2]/3.),
        np.add(np.add(m[0]/1.5, m[1]/1.5), m[2]/3.),
        np.add(np.add(m[0]/3., m[1]/3.), m[2]/1.5),
        np.add(np.add(m[0]/1.5, m[1]/3.), m[2]/1.5),
        np.add(np.add(m[0]/3., m[1]/1.5), m[2]/1.5),
        np.add(np.add(m[0]/1.5, m[1]/1.5), m[2]/1.5)
    ])
    cell_center = np.mean(cell_vertices, axis=0)
    other_cell_centers = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                c = np.add(cell_center, np.multiply(i, m_0[0]))
                c = np.add(c, np.multiply(j, m_0[1]))
                c = np.add(c, np.multiply(k, m_0[2]))
                other_cell_centers.append(c)

    max_distance_in_cell = sq_dist(cell_vertices[0], cell_center)

    points = [s.coords for s in st.sites]
    radii = [float(s.specie.atomic_radius) for s in st.sites]

    # Create the initial Delaunay triangulation of all sites in the
    # supercell.
    delaunay = Delaunay(points)
    all_simplices = delaunay.simplices.copy()

    # Now filter those Delaunay simplices to only those with
    # at least one vertex lying within the center unit cell.
    simplices = []
    center_cell = ConvexHull(cell_vertices)
    if not octahedra:
        for simplex in all_simplices:
            for vertex in simplex:
                if sq_dist(cell_center, points[vertex]) <= max_distance_in_cell\
                        and sq_dist(cell_center, points[vertex]) ==\
                        min([sq_dist(points[vertex], pt) for pt in
                             other_cell_centers]):
                    simplices.append(simplex)
                    break
    else:
        for simplex in all_simplices:
            n = 0
            for vertex in simplex:
                if sq_dist(cell_center, points[vertex]) <= max_distance_in_cell\
                        and sq_dist(cell_center, points[vertex]) ==\
                        min([sq_dist(points[vertex], pt) for pt in
                             other_cell_centers]):
                    n += 1
            if n == 4:
                simplices.append(simplex)

    # Calculate the maximum interstitial
    # radius for all the relevant tetrahedra.
    tetrahedra = []
    for simplex in simplices:
        a = points[simplex[0]]
        r_a = radii[simplex[0]]
        b = points[simplex[1]]
        r_b = radii[simplex[1]]
        c = points[simplex[2]]
        r_c = radii[simplex[2]]
        d = points[simplex[3]]
        r_d = radii[simplex[3]]
        centroid = np.mean([a,b,c,d], axis=0)

        # Add the atomic radii to the nuclei loactions to find
        # their "true" extrema, then use these to find the
        # "true" centroid.
        move = 1
        while move > 0.01:
            true_a = pt_btwn(a, centroid, r_a)
            true_b = pt_btwn(b, centroid, r_b)
            true_c = pt_btwn(c, centroid, r_c)
            true_d = pt_btwn(d, centroid, r_d)
            true_centroid = np.mean(
                [true_a,true_b,true_c,true_d], axis=0
            )
            move = sq_dist(true_centroid, centroid)
            centroid = true_centroid

        max_radius = sqrt(min(
            [sq_dist(true_centroid, pt) for pt in [true_a,true_b,true_c,true_d]]
        ))

        tetrahedra.append(
            (true_centroid, [tuple(x) for x in [a, b, c, d]],
             [r_a, r_b, r_c, r_d], 4, max_radius)
        )

    interstitials = {"tetrahedral": []}
    if octahedra:
        tet_pts = [i[1] for i in tetrahedra]
        tet_pts = list(set([coords for pt in tet_pts for coords in pt]))
        interstitials.update({"hexahedral": [], "octahedral": []})
        for i in range(len(tetrahedra)):
            for j in range(i, len(tetrahedra)):
                # If 3 vertices are shared then the tetrahedra
                # share a face and form a bipyramid.
                shared = list(set(tetrahedra[i][1]) & set(tetrahedra[j][1]))
                if len(shared) == 3:
                    # Vertices of the bipyramid
                    a = tetrahedra[i][1][0]
                    r_a = tetrahedra[i][2][0]
                    b = tetrahedra[i][1][1]
                    r_b = tetrahedra[i][2][1]
                    c = tetrahedra[i][1][2]
                    r_c = tetrahedra[i][2][2]
                    d = tetrahedra[i][1][3]
                    r_d = tetrahedra[i][2][3]
                    # Fifth point to define trigonal bipyramid
                    e, r_e = [
                        (s, tetrahedra[j][2][k]) for k, s in
                        enumerate(tetrahedra[j][1]) if s
                        not in tetrahedra[i][1]
                    ][0]

                    h_centroid = np.mean([a, b, c, d, e], axis=0)
                    move = 1
                    while move > 0.01:
                        true_a = pt_btwn(a, h_centroid, r_a)
                        true_b = pt_btwn(b, h_centroid, r_b)
                        true_c = pt_btwn(c, h_centroid, r_c)
                        true_d = pt_btwn(d, h_centroid, r_d)
                        true_e = pt_btwn(e, h_centroid, r_e)

                        true_h_centroid = np.mean(
                            [true_a,true_b,true_c,true_d,true_e], axis=0
                        )
                        move = sq_dist(true_h_centroid, h_centroid)
                        h_centroid = true_h_centroid

                    r_h = sqrt(min(
                        [sq_dist(true_h_centroid, pt) for pt in
                         [true_a, true_b, true_c, true_d, true_e]]
                    ))

                    # Add the bipyramid to the final list
                    # of interstitials.
                    interstitials["hexahedral"].append(
                        (tuple(h_centroid), r_h)
                    )

                    # Enlarge the bipyramid by one point to create
                    # octahedra.
                    v1 = np.subtract(shared[0], shared[1])
                    v2 = np.subtract(shared[0], shared[2])
                    tol = max([sq_dist(shared[0], shared[1]),
                               sq_dist(shared[0], shared[2]),
                               sq_dist(shared[1], shared[2])]) * 1.1
                    for index, f in enumerate(tet_pts):
                        v3 = np.subtract(shared[0], f)
                        distances = [sq_dist(f, p) for p in shared]
                        distances.sort()
                        if 0 < distances[0] < tol and 0 < distances[1] < tol\
                                and np.dot(v3, (np.cross(v1, v2))) == 0:
                            r_f = radii[index]
                            o_centroid = np.mean([a, b, c, d, e, f], axis=0)

                            move = 1
                            while move > 0.01:
                                true_a = pt_btwn(a, o_centroid, r_a)
                                true_b = pt_btwn(b, o_centroid, r_b)
                                true_c = pt_btwn(c, o_centroid, r_c)
                                true_d = pt_btwn(d, o_centroid, r_d)
                                true_e = pt_btwn(e, o_centroid, r_e)
                                true_f = pt_btwn(f, o_centroid, r_f)

                                true_o_centroid = np.mean(
                                    [true_a,true_b,true_c,true_d,true_e,true_f],
                                    axis=0
                                )
                                move = sq_dist(true_o_centroid, o_centroid)
                                o_centroid = true_o_centroid

                            r_o = sqrt(min(
                                [sq_dist(true_o_centroid, pt) for
                                 pt in [true_a,true_b,true_c,true_d,true_e,
                                        true_f]]
                            ))

                            # Add the octahedron to the final
                            # list of interstitials.
                            interstitials["octahedral"].append(
                                (tuple(o_centroid), r_o)
                            )
        interstitials["hexahedral"] = list(set(interstitials["hexahedral"]))
        interstitials["octahedral"] = list(set(interstitials["octahedral"]))


    interstitials["tetrahedral"] = [(i[0], i[4]) for i in tetrahedra]

    # Since the centroid coordinates were given in the center
    # cell of the supercell, bring them back into the original
    # unit cell.
    if n_sites < 4:
        f = 1./3.
    else:
        f = 1.
    for c in interstitials:
        for i in range(len(interstitials[c])):
            for r in m_0:
                interstitials[c][i] = (
                    np.multiply(
                        np.subtract(np.array(interstitials[c][i][0]), r), f
                    ),
                    interstitials[c][i][1]
                )

    # Sort by the maximum radii
    for c in interstitials:
        interstitials[c].sort(key=operator.itemgetter(1))
        interstitials[c].reverse()

    if unique:
        sga = SpacegroupAnalyzer(structure)
        sop = sga.get_space_group_operations()
        l = structure.lattice
        for c in interstitials:
            remove = []
            for i in range(len(interstitials[c])):
                if i not in remove:
                    site_i = PeriodicSite("C", interstitials[c][i][0], l)
                    for j in range(i+1, len(interstitials[c])):
                        if interstitials[c][i][1] == interstitials[c][j][1] and\
                                sop.are_symmetrically_equivalent(
                                    [site_i],
                                    [PeriodicSite("C",interstitials[c][j][0],l)]
                                ):
                            remove.append(j)
            interstitials[c] = [interstitials[c][x] for x in
                                range(len(interstitials[c])) if x not in remove]

    return interstitials


def get_coordination_polyhedra(structure, cation, anion="O"):

    r_c, r_a = Element(cation).atomic_radius, Element(anion).atomic_radius

    st = structure.copy()
    cations = [s for s in st.sites if s.specie.symbol == cation]
    uc_tetrahedra, uc_octahedra = [], []
    for s in cations:
        anion_shell = [a[0] for a in st.get_neighbors(s, (r_c+r_a)*1.1)]
        if len(anion_shell) == 4:
            uc_tetrahedra.append(
                [tuple([round(c, 3) for c in a.coords]) for a in anion_shell])
        elif len(anion_shell) == 6:
            uc_octahedra.append(
                [tuple([round(c, 3) for c in a.coords]) for a in anion_shell])

    st.make_supercell(2)
    cations = [s for s in st.sites if s.specie.symbol == cation]
    tetrahedra, octahedra = [], []
    for s in cations:
        anion_shell = [a[0] for a in st.get_neighbors(s, (r_c+r_a)*1.1)]
        if len(anion_shell) == 4:
            tetrahedra.append(
                [tuple([round(c, 3) for c in a.coords]) for a in anion_shell])
        elif len(anion_shell) == 6:
            octahedra.append(
                [tuple([round(c, 3) for c in a.coords]) for a in anion_shell])

    t_corner, t_edge, t_face = [], [], []
    o_corner, o_edge, o_face = [], [], []
    if len(tetrahedra) != 0:
        for i in range(len(tetrahedra)):
            t1 = tetrahedra[i]
            for j in range(i+1, len(tetrahedra)):
                t2 = tetrahedra[j]
                shared = list(set(t1) & set(t2))
                if len(shared) == 1:
                    # Corner sharing
                    if t1 in uc_tetrahedra and t1 not in t_corner:
                        t_corner.append(t1)
                    if t2 in uc_tetrahedra and t2 not in t_corner:
                        t_corner.append(t2)
                elif len(shared) == 2:
                    # Edge sharing
                    if t1 in uc_tetrahedra and t1 not in t_edge:
                        t_edge.append(t1)
                    if t2 in uc_tetrahedra and t2 not in t_edge:
                        t_edge.append(t2)
                elif len(shared) == 3:
                    # Face sharing
                    if t1 in uc_tetrahedra and t1 not in t_face:
                        t_face.append(t1)
                    if t2 in uc_tetrahedra and t2 not in t_face:
                        t_face.append(t2)
    if len(octahedra) != 0:
        for i in range(len(octahedra)):
            o1 = octahedra[i]
            for j in range(i+1, len(octahedra)):
                o2 = octahedra[j]
                shared = list(set(o1) & set(o2))
                if len(shared) == 1:
                    # Corner sharing
                    if o1 in uc_octahedra and o1 not in o_corner:
                        o_corner.append(o1)
                    if o2 in uc_octahedra and o2 not in o_corner:
                        o_corner.append(o2)
                elif len(shared) == 2:
                    # Edge sharing
                    if o1 in uc_octahedra and o1 not in o_edge:
                        o_edge.append(o1)
                    if o2 in uc_octahedra and o2 not in o_edge:
                        o_edge.append(o2)
                elif len(shared) == 3:
                    # Face sharing
                    if o1 in uc_octahedra and o1 not in o_face:
                        o_face.append(o1)
                    if o2 in uc_octahedra and o2 not in o_face:
                        o_face.append(o2)

    polyhedra = {
        "tetrahedra": {"corner": t_edge, "edge": t_corner, "face": t_face},
        "octahedra": {"corner": o_edge, "edge": o_corner, "face": o_face}
    }
    return polyhedra



def plot_ion_hull_and_voltages(ion, charge=None, fmt='pdf'):
    """
    Plots the phase diagram between the pure material and pure ion,
    Connecting the points on the convex hull of the phase diagram.
    Args:
        ion (str): name of atom that was intercalated, e.g. 'Li'.
        charge (float): charge donated by each ion.
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    Returns:
        capacity (float): Maximum capacity
    """

    # Calculated with the relax() function in
    # mat2d.stability.startup. If you are using other input
    # parameters, you need to recalculate these values!
    ion_ev_fu = {'Li': -1.838, 'Mg': 0.620, 'Al': -3.291}

    if charge is None:
        charge = Element(ion).common_oxidation_states[0]

    energy = Vasprun('vasprun.xml').final_energy
    composition = Structure.from_file('POSCAR').composition

    # Get the formula (with single-digit integers preceded by a '_').
    twod_material = list(composition.reduced_formula)
    twod_formula = str()
    for i in range(len(twod_material)):
        try:
            int(twod_material[i])
            twod_formula += '_{}'.format(twod_material[i])
        except:
            twod_formula += twod_material[i]

    twod_ev_fu = energy / composition.get_reduced_composition_and_factor()[1]

    data = [(0, 0, 0, twod_ev_fu)]  # (at% ion, n_ions, E_F, abs_energy)
    dirs = [dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir)]
    for directory in dirs:
        if is_converged(directory):
            os.chdir(directory)
            energy = Vasprun('vasprun.xml').final_energy
            composition = Structure.from_file('POSCAR').composition
            ion_fraction = composition.get_atomic_fraction(ion)

            no_ion_comp_dict = composition.as_dict()
            no_ion_comp_dict.update({ion: 0})
            no_ion_comp = Composition.from_dict(no_ion_comp_dict)

            n_twod_fu = no_ion_comp.get_reduced_composition_and_factor()[1]
            n_ions = composition[ion] / n_twod_fu

            E_F = ((energy - composition[ion] * ion_ev_fu[ion] -
                    twod_ev_fu * n_twod_fu)/ composition.num_atoms)

            data.append((ion_fraction, n_ions, E_F, energy / n_twod_fu))

            os.chdir('../')
    data.append((1, 1, 0, ion_ev_fu[ion]))  # Pure ion

    sorted_data = sorted(data, key=operator.itemgetter(0))

    # Determine which compositions are on the convex hull.
    energy_profile = np.array([[item[0], item[2]] for item in sorted_data
                               if item[2] <= 0])
    hull = ConvexHull(energy_profile)
    convex_ion_fractions = [energy_profile[vertex, 0] for vertex in hull.vertices]
    convex_formation_energies = [energy_profile[vertex, 1] for vertex in hull.vertices]

    convex_ion_fractions.append(convex_ion_fractions.pop(0))
    convex_formation_energies.append(convex_formation_energies.pop(0))

    concave_ion_fractions = [pt[0] for pt in sorted_data
                             if pt[0] not in convex_ion_fractions]
    concave_formation_energies = [pt[2] for pt in sorted_data
                                  if pt[0] not in convex_ion_fractions]

    for item in data:
        if item[0] == sorted(convex_ion_fractions)[-2]:
            max_ions = item[1]
            molar_mass = Composition(no_ion_comp.reduced_formula).weight
            faraday = 26801  # In mAh/mol
            capacity = (max_ions * charge * faraday) / molar_mass  # In mAh/g

    voltage_profile = []
    j = 0
    k = 0
    for i in range(1, len(sorted_data) - 1):
        if sorted_data[i][0] in convex_ion_fractions:
            voltage = -(((sorted_data[i][3] - sorted_data[k][3])-
                         (sorted_data[i][1] - sorted_data[k][1]) * ion_ev_fu[ion])
                        / (sorted_data[i][1] - sorted_data[k][1]))
            voltage_profile.append((sorted_data[k][0], voltage))
            voltage_profile.append((sorted_data[i][0], voltage))
            j += 1
            k = i

    voltage_profile.append((voltage_profile[-1][0], 0))
    voltage_profile.append((1, 0))

    voltage_profile_x = [tup[0] for tup in voltage_profile]
    voltage_profile_y = [tup[1] for tup in voltage_profile]

    ax = plt.figure(figsize=(14, 10)).gca()

    ax.plot([0, 1], [0, 0], 'k--')
    ax.plot(convex_ion_fractions, convex_formation_energies, 'b-', marker='o',
            markersize=12, markeredgecolor='none')
    ax.plot(concave_ion_fractions, concave_formation_energies, 'r', marker='o',
            linewidth=0, markersize=12, markeredgecolor='none')

    ax2 = ax.twinx()
    ax2.plot(voltage_profile_x, voltage_profile_y, 'k-', marker='o')

    ax.text(0, 0.002, r'$\mathrm{%s}$' % twod_formula, family='serif', size=24)
    ax.text(0.99, 0.002, r'$\mathrm{%s}$' % ion, family='serif', size=24,
            horizontalalignment='right')

    ax.set_xticklabels(ax.get_xticks(), family='serif', size=20)
    ax.set_yticklabels(ax.get_yticks(), family='serif', size=20)
    ax2.set_yticklabels(ax2.get_yticks(), family='serif', size=20)

    ax.set_xlabel('at% {}'.format(ion), family='serif', size=28)
    ax.set_ylabel(r'$\mathrm{E_F\/(eV/atom)}$', size=28)

    ax2.yaxis.set_label_position('right')
    if ion == 'Li':
        ax2.set_ylabel(r'$\mathrm{Potential\/vs.\/Li/Li^+\/(V)}$', size=28)
    elif ion == 'Mg':
        ax2.set_ylabel(r'$\mathrm{Potential\/vs.\/Mg/Mg^{2+}\/(V)}$', size=28)
    elif ion == 'Al':
        ax2.set_ylabel(r'$\mathrm{Potential\/vs.\/Al/Al^{3+}\/(V)}$', size=28)

    plt.savefig('{}_hull.{}'.format(ion, fmt), transparent=True)

    return capacity  # In mAh/g
