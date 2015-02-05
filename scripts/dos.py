from __future__ import division, unicode_literals, print_function

"""
reads in vasprun.xml file and plots the density of states 
"""
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

if __name__ == "__main__":
    # readin the density of states from vasprun.xml file
    run = Vasprun("vasprun.xml", parse_projected_eigen = True)
    complete_dos = run.complete_dos
    print('cbm and vbm ', complete_dos.get_cbm_vbm())
    print('gap = ', complete_dos.get_gap())
    # get orbital projected DOS.    
    spd_dos = complete_dos.get_spd_dos()
    plotter = DosPlotter()
    plotter.add_dos_dict(spd_dos)
    plotter.save_plot('dos.eps')
