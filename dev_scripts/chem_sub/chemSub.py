#Written by Joshua Paul
#E-mail: joshua.thomas.paul@gmail.com
#Last Updated: 1/24/2017

#   This script performs chemical substitutions for binary systems. It will take any number of monolayers
#   and their corresponding bulk structures (required) and perform chemical substitution on them using two
#   lists of elements (Ms and Xs, below). It will also find the competing phase(s) for the monolayer,
#   whether that is a structure with the same ratio of elements as the monolayers (aka

#   In order to run this script, you will need to create a directory called 'inputs' in your home directory.
#   Inside of it, you will place the INCAR you wish to use for all of your simulations, the submit file for
#   2D materials (submit_2D) and the submit file for 3D materials (submit). Next you will need to create
#   two directories: monolayers, and bulks. Inside monolayers you will place the monolayer structures you
#   wish to simulate, vacuum padding not required. However, the vacuum will be applied along the c-axis.
#   Inside bulks, you will place the bulk structure that each monolayer came from, with the same names.
#   For example, if you put a POSCAR in monolayers and name it '1T', the bulk structure that the monolayer
#   came from will also need to be named '1T'.
#   To run bader charge analysis, you must place the 'bader' and 'chgsum.pl' scripts in your bin directory.

import os

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.phasediagram.maker import PhaseDiagram

from mpinterfaces import MPR, get_struct_from_mp

__author__ = "Joshua T. Paul, Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

#home = os.path.expanduser('~')
def get_competing_phases_new(structure):
    """
    Collect the species to which the 2D materials might decompose to.
    Since a lot of 2D materials with similar compositions will have the
    same competing phases, duplicates aren't counted.
    """

    total_competing_phases = []
    composition = structure.composition
    energy = 100
    my_entry = ComputedEntry(composition, energy)  # 2D material
    entries = MPR.get_entries_in_chemsys(elements=[elt.symbol for elt in composition])

    entries.append(my_entry)  # 2D material

    pda = PDAnalyzer(PhaseDiagram(entries))
    decomp = pda.get_decomp_and_e_above_hull(my_entry, allow_negative=True)
    competing_phases = [
        (entry.composition.reduced_formula,
         entry.entry_id) for entry in decomp[0]
        ]

    # Keep a running list of all unique competing phases, since in
    # high throughput 2D searches there is usually some overlap in
    # competing phases for different materials.
    for specie in competing_phases:
        if specie not in total_competing_phases:
            total_competing_phases.append(specie)

    return total_competing_phases


def makePoscars(M_elements, X_elements, monolayer_path, bulk_path, getCompeting = True):
    """
    function that writes out poscars of 
    chem-substituted M-X monolayers

    Args:
       M_elements : (list) M
       X_elements : (list) X
       ratio      : (float?)  
       monolayers : (str) path to poscar files of monolayers
       bulks      : (str) path to poscar files of bulk 

    Returns:
       List of Poscar file objects for twod monolayers and their competitors 
    """
    # checing the existence of the monolayer
    # and bulk directories 
    # these direcotries can simply contain the motifs at 
    # the top most level 
    monosDir = os.listdir(monolayer_path)
    bulksDir = os.listdir(bulk_path)
    # if bulks in os.listdir('inputs'):
    #    bulksDir = os.listdir('inputs/bulks')
    # else:
    #    bulks = []
    twod = [] 
    competitors = [] 
    monos = {}
    bulks = {}

    results = {}

    # making direcotries of types of chem sub directories 
    results['onRatio'] = []
    results['hull'] = []
    # read in monolayer and bulk motifs
    for mon in monosDir:
        monos[mon] = Structure.from_file(monolayer_path+'/'+mon)
        results[mon] = []
    for bulk in bulksDir:
        bulks[bulk] = Structure.from_file(bulk_path+'/'+bulk)
        results[bulk+'_bulk'] = []
    # perform the chem sub based on the input X and M lists 
    for base_name, mono in monos.items():
        # os.mkdir(monolayer)
        # os.chdir(monolayer)
        for M in M_elements:
            for X in X_elements:
                coreElements = [M,X]
                # naming the files .. can be set as poscar.comment
                sub_spec = np.unique(mono.species)
                if len(sub_spec) != len(coreElements):
                    print (base_name+' DOES NOT HAVE '+len(coreElements)+' ELEMENTS. ENDING LOOP')
                    break

                else:
                    substitutionDict = {}
                    counter = 0
                    for element in sub_spec:
                        substitutionDict[element] = coreElements[counter]
                        counter +=1
                    substitutedStruct = mono
                    substitutedStruct.replace_species(substitutionDict)                    
                    results[base_name].append(substitutedStruct)
               
                

                if getCompeting:
                    competing = get_competing_phases_new(mono) # argument taken is the POSCAR in the current directory.
                                                       # Modify get_competing_phases code to allow input structure 
                                                       # object?
                    if len(competing) == 1:
                        # ONRATIO
                        newStruct = get_struct_from_mp(competing[0][1]) # assume the  
                                                                         # lowest hull distance 
                        
                        results['onRatio'].append(newStruct)
                                                         # chem sub 
                   ## how about more than 2 ? this would just be the case in ternaries, etc.
                    elif len(competing) == 2:
                        # HULL 
                        newStruct1 = get_struct_from_mp(competing[0][1])
                        newStruct2 = get_struct_from_mp(competing[1][1])
                        results['hull'].append(newStruct1)
                        results['hull'].append(newStruct2)
                        # write_potcar(pathToPOTCAR)
                        # relax(submit=False)
    # same for bulks 
    for base_name, bulk in bulks.items():
        for M in M_elements:
            for X in X_elements:
                coreElements = [M,X]
                sub_spec = np.unique(bulk.species)
                if len(sub_spec) != len(coreElements):
                    print (base_name+' DOES NOT HAVE '+len(coreElements)+' ELEMENTS. ENDING LOOP')
                    break
                else:
                    substitutionDict = {}
                    counter = 0
                    for element in sub_spec:
                        substitutionDict[element] = coreElements[counter]
                        counter+=1
                    substitutedStruct = bulk
                    substitutedStruct.replace_species(substitutionDict)
                    results[base_name+'_bulk'].append(substitutedStruct)


    return results
if __name__=='__main__':
    Ms = ['Ir']
    Xs = ['I','Cl']
    print (makePoscars(Ms, Xs, 'monos','bulks'))




