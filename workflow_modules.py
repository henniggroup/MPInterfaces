import os

import vasp_tools as vt

from pymatgen.io.vasp.inputs import Incar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.matproj.rest import MPRester

from maker_2d import Pourbaix2D

from matplotlib import cm
import matplotlib.pyplot as plt

from monty.serialization import loadfn


INCAR_DICT = {}
HSE_INCAR_DICT = {}
KERNEL_PATH = '/home/mashton/vdw_kernel.bindat'
MPR = MPRester(
    loadfn(os.path.join(os.path.expanduser('~', 'dbauth.yaml')))['mp_api']
    )


def relax(directories, submit=True):

    for directory in directories:
        os.chdir(directory)

        # Ensure 20A interlayer vacuum
        vt.add_vacuum(20 - vt.get_spacing(), 0.9)

        # vdw_kernel.bindat file required for VDW calculations.
        os.system('cp /home/mashton/vdw_kernel.bindat .')

        # KPOINTS
        Kpoints.automatic_density(Structure.from_file('POSCAR'),
                                  1000).write_file('KPOINTS')

        # INCAR
        Incar.from_dict(INCAR_DICT).write_file('INCAR')

        # POTCAR
        vt.write_potcar()

        # Submission script
        vt.write_runjob(directory, 1, 8, '600mb', '6:00:00', 'vasp_noz')

        if submit:
            os.system('qsub runjob')
        os.chdir('../')


def is_converged(directory):

    if 'reached required accuracy' in open(
        '{}/job.log'.format(directory)
            ).read():
        return True
    else:
        return False


def get_competing_species(directories):

    total_competing_species = []
    hull_distances = {}
    finished_competitors = {}

    # Determine which competing species have been relaxed in the current
    # framework and store them in a dictionary ({formula: entry}).
    if os.path.isdir('all_competitors'):
        os.chdir('all_competitors')
        for comp_dir in [
            dir for dir in os.listdirs(os.getcwd()) if os.path.isdir(dir) and
            is_converged(dir)
                ]:
            os.chdir(comp_dir)
            composition = Structure.from_file('POSCAR').composition
            energy = Vasprun('vasprun.xml').final_energy
            finished_competitors[comp_dir] = ComputedEntry(composition, energy)
            os.chdir('../')
        os.chdir('../')

    for directory in directories:
        os.chdir(directory)
        composition = Structure.from_file('POSCAR').composition
        energy = Vasprun('vasprun.xml').final_energy
        my_entry = ComputedEntry(composition, energy)  # 2D material
        entries = MPR.get_entries_in_chemsys(
            [elt.symbol for elt in composition]
            )

        # If the energies of competing species have been calculated in
        # the current framework, put them in the phase diagram instead
        # of the MP energies.
        for i in range(len(entries)):
            formula = entries[i].composition.reduced_formula
            if formula in finished_competitors:
                entries[i] = finished_competitors[formula]

        entries.append(my_entry)  # 2D material

        pda = PDAnalyzer(PhaseDiagram(entries))
        decomp = pda.get_decomposition_and_e_above_hull(my_entry,
                                                        allow_negative=True)
        competing_species = [
            (entry.composition.reduced_formula,
             entry.entry_id) for entry in decomp[0]
            ]

        # Keep a running list of all unique competing species, since in
        # high throughput 2D searches there is usually some overlap in
        # competing species for different materials.
        for specie in competing_species:
            if specie not in total_competing_species:
                total_competing_species.append(specie)
        os.chdir('../')
        hull_distances[composition.reduced_formula] = decomp[1]

    return [total_competing_species, hull_distances]


def relax_competing_species(competing_species, submit=True):

    if not os.path.isdir('all_competitors'):
        os.mkdir('all_competitors')
    os.chdir('all_competitors')

    for specie in competing_species:
        if not os.path.isdir(specie):
            os.mkdir(specie[0])
        os.chdir(specie[0])
        os.system('cp /home/mashton/vdw_kernel.bindat .')
        structure = MPR.get_structure_by_material_id(specie[1])
        structure.to('POSCAR', 'POSCAR')
        Kpoints.automatic_density(structure, 1000).write_file('KPOINTS')
        Incar.from_dict(INCAR_DICT).write_file('INCAR')
        vt.write_potcar()
        vt.write_runjob(specie[0], 1, 8, '600mb', '6:00:00', 'vasp')
        if submit:
            os.system('qsub runjob')
        os.chdir('../')
    os.chdir('../')


def is_pourbaix_stable(directory):

    composition = Structure.from_file('POSCAR').composition
    energy = Vasprun('vasprun.xml').final_energy
    # Pourbaix2D requires the entry to be in a database right now- that
    # should change.
    pourbaix_plotter = Pourbaix2D(composition, energy, metastability=0)
    # This function doesn't exist yet- I need to write it.
    return pourbaix_plotter.is_stable()


def get_band_structures(directories):

    band_gaps = {}
    for directory in directories:
        os.chdir(directory)
        vasprun = Vasprun('vasprun.xml')
        band_gap = vasprun.get_band_structure().get_band_gap()

        transition = band_gap['transition'].split('-')

        if transition[0] == transition[1]:
            is_direct = True
        else:
            is_direct = False

        if band_gap['energy']:
            cbm = vasprun.get_band_structure().get_cbm()
            vbm = vasprun.get_band_structure().get_vbm()
            band_gaps[directory] = {'CBM': cbm, 'VBM': vbm,
                                    'Direct': is_direct}
        else:
            band_gaps[directory] = False

    return band_gaps


def run_hse_calculations(directories, submit=True):

    for directory in directories:
        os.chdir(directory)
        vasprun = Vasprun('vasprun.xml')
        band_gap = vasprun.get_band_structure().get_band_gap()
        kpath = []
        if band_gap['energy']:
            transition = band_gap['transition'].split('-')

            increment = ((transition[1][0] - transition[0][0]) / 9,
                         (transition[1][1] - transition[0][1]) / 9,
                         (transition[1][2] - transition[0][2] / 9))

            for i in range(10):
                kpath.append((transition[0] + increment[0] * i),
                             (transition[1] + increment[1] * i),
                             (transition[2] + increment[2] * i))

            if not os.path.isdir('HSE'):
                os.mkdir('HSE')
            os.chdir('HSE')
            os.system('cp ../CONTCAR ./POSCAR')
            os.system('cp ../POTCAR ./POTCAR')
            os.system('cp ../vdw_kernel.bindat ./')
            Incar.from_dict(HSE_INCAR_DICT).write_file('INCAR')
            vt.write_runjob(directory, 1, 32, '1200mb', '150:00:00', 'vasp')

            os.system('cp ../IBZKPT ./KPOINTS')
            kpoints_lines = open('KPOINTS').readlines()
            n_kpts = int(kpoints_lines[1].split()[0])
            with open('KPOINTS', 'w') as kpts:
                kpts.write(kpoints_lines[0])
                kpts.write('{}\n'.format(n_kpts + 10))
                for line in kpoints_lines[2:]:
                    kpts.write(line)
                for kpoint in kpath:
                    kpts.write('{}\n'.format(' '.join(kpoint)))

            if submit:
                os.system('qsub runjob')
        os.chdir('../')


def run_friction_calculations(directories, submit=True):

    for directory in directories:
        os.chdir(directory)
        vt.add_vacuum(3 - vt.get_spacing(), 0.8)
        structure = Structure.from_file('POSCAR')
        n_sites_per_layer = structure.num_sites
        structure.make_supercell([1, 1, 2])
        structure.to('POSCAR', 'POSCAR')
        vt.add_vacuum(12, 0.9)

        z_coords = []
        for site in structure.sites:
            z_coords.append(site.z)
        bottom_layer_max_height = sorted(z_coords)[n_sites_per_layer - 1]

        for x in range(10):
            for y in range(10):
                dir = '{}x{}'.format(x, y)

                # Copy input files
                os.system('cp INCAR {}/'.format(dir))
                os.system('cp KPOINTS {}/'.format(dir))
                os.system('cp POSCAR {}/'.format(dir))
                os.system('cp POTCAR {}/'.format(dir))
                os.system('cp vdw_kernel.bindat {}/'.format(dir))

                os.chdir(dir)
                incar_dict = Incar.from_file('INCAR').as_dict()
                incar_dict.update({'NSW': 0})
                Incar.from_dict(incar_dict).write_file('INCAR')

                # Shift the top layer
                poscar_lines = open('POSCAR').readlines()
                with open('POSCAR', 'w') as poscar:
                    for line in poscar_lines[:8]:
                        poscar.write(line)
                    for line in poscar_lines[8:8 + structure.num_sites]:
                        split_line = line.split()
                        if float(split_line[2]) > bottom_layer_max_height:
                            new_coords = [float(split_line[0]) + float(x)/10.0,
                                          float(split_line[1]) + float(y)/10.0,
                                          float(split_line[2])]
                            poscar.write(' '.join([str(i) for i in new_coords])
                                         + '\n')
                        else:
                            poscar.write(line)

                vt.write_runjob(dir, 1, 8, '400mb', '1:00:00', 'vasp')

                if submit:
                    os.system('qsub runjob')

                os.chdir('../')
        os.chdir('../')


def plot_hull_distances(hull_distances):

    fig = plt.figure()
    ax = fig.gca()
    ax.set_ylim(0, 1000)
    ax.set_xlim(0, len(hull_distances))

    x_ticklabels = []
    i = 0
    for compound in hull_distances:
        x_ticklabels.append(compound)
        if hull_distances[compound] < 100:
            color_code = 0.5
        elif hull_distances[compound] < 200:
            color_code = 0.71
        else:
            color_code = 0.92

        ax.add_patch(plt.Rectangle(i, 0), height=hull_distances[compound],
                     width=1, linewidth=0, facecolor=cm.jet(color_code))

    ax.set_xticklabels(x_ticklabels)

    plt.savefig('stability_plot.pdf', transparent=True)


def plot_friction_surface(directory):

    ax = plt.figure(figsize=(10, 10)).gca()

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)

    ENERGY_ARRAY = []

    X_VALUES = range(10)
    Y_VALUES = range(10)

    for x in X_VALUES:
        ENERGY_ARRAY.append([])
        for y in Y_VALUES:
            dir = '{}x{}'.format(x, y)
            os.chdir(dir)
            ENERGY_ARRAY[x].append(vt.get_toten())
            os.chdir('../')

    minima = []
    maxima = []
    for x in X_VALUES:
        minima.append(min(ENERGY_ARRAY[x]))
        maxima.append(max(ENERGY_ARRAY[x]))
    abs_minimum = min(minima)
    abs_maximum = max(maxima)
    print abs_minimum, abs_maximum
    print ENERGY_ARRAY[9][9]

    for x in X_VALUES:
        for y in Y_VALUES:
            scaled_energy = ENERGY_ARRAY[x][y] - abs_minimum
            color_code = scaled_energy / (abs_maximum - abs_minimum)

            ax.add_patch(plt.Rectangle((x, y), width=1, height=1,
                                       facecolor=cm.coolwarm(color_code),
                                       linewidth=0))

    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_ticks([])

    plt.savefig('{}.pdf'.format(os.getcwd().split('/')[-1]),
                transparent='True')


def plot_band_alignments_relative_to_CO2(band_gaps):

    ax = plt.figure(figsize=(10, 10)).gca()
    ax.set_xlim(len(band_gaps))

    x_ticklabels = []
    vbms = []
    i = 0
    for compound in band_gaps:
        x_ticklabels.append(compound)
        cbm = band_gaps[compound]['CBM']
        vbm = band_gaps[compound]['VBM']
        vbms.append(vbm)
        if band_gaps[compound]['Direct']:
            linewidth = 3
        else:
            linewidth = 0

        ax.add_patch(plt.Rectangle((i, cbm), height=-cbm, width=0.9,
                                   facecolor="#0066cc", linewidth=linewidth))
        i += 1

    ax.set_ylim(min(vbms)-0.5, 0)

    ax.plot([0, i], [-1.9, -1.9], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -1.94, r'$\mathrm{CO_2+\/e^-\/\rightarrow\/CO^-_2}$',
            size=20)

    ax.plot([0, i], [-0.61, -0.61], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.7,
            r'$\mathrm{CO_2\/+\/2H^+\/+\/2e^-\/\rightarrow\/HCO_2H}$', size=20)

    ax.plot([0, i], [-0.53, -0.53], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.59,
            r'$\mathrm{CO_2\/+\/2H^+\/+\/2e^-\/\rightarrow\/CO\/+\/H_2O}$',
            size=20)

    ax.plot([0, i], [-0.48, -0.48], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.48,
            r'$\mathrm{CO_2\/+\/4H^+\/+\/4e^-\/\rightarrow\/HCHO\/+\/H_2O}$',
            size=20)

    ax.plot([0, i], [-0.38, -0.38], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.37,
            r'$\mathrm{CO_2\/+\/6H^+\/+\/6e^-\/\rightarrow\/CH_3OH\/+\/H_2O}$',
            size=20)

    ax.plot([0, 16], [-0.24, -0.24], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.24,
            r'$\mathrm{CO_2\/+\/8H^+\/+\/8e^-\/\rightarrow\/CH_4\/+\/2H_2O}$',
            size=20)

    ax.set_xticks(range(i))
    ax.set_xticklabels(x_ticklabels, family='serif', size=20, rotation=60)
    ax.set_yticklabels(ax.get_yticks(), family='serif', size=20)

    ax.add_patch(plt.Rectangle((i + 1, min(vbms) + 0.25), width=3, height=0.2,
                               facecolor='w', linewidth=3))
    ax.text(i + 2.5, min(vbms) + 0.35, 'Direct', family='serif', color='k',
            size=20, horizontalalignment='center', verticalalignment='center')
    ax.add_patch(plt.Rectangle((i + 5, min(vbms) + 0.25), width=3, height=0.2,
                               facecolor='w', linewidth=0))
    ax.text(i + 6.5, min(vbms) + 0.35, 'Indirect', family='serif', size=20,
            color='k', horizontalalignment='center',
            verticalalignment='center')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.savefig('band_alignments.pdf', transparent=True)
