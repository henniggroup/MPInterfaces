import os

from twod_materials.utils import write_runjob

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Incar


HSE_INCAR_DICT = {}


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
            write_runjob(directory, 1, 32, '1200mb', '150:00:00', 'vasp')

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
