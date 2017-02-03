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
import sys
import time
import collections
import shutil
from pymatgen.io.vasp.iputs import Poscar
from pymatgen.core.structure import Structure
from mpinterfaces import get_struct_from_mp
from mpinterfaces.utils import add_vacuum_padding
from twod_materials.utils import write_potcar # test of this is redundant or there is an extra feature 
# from twod_materials.stability.startup import relax
from twod_materials.stability.analysis import get_competing_phases

import numpy

#home = os.path.expanduser('~')

def makePoscars(M_elements, X_elements, ratio, mag=False, monolayer_path, bulk_path):
    nameOfErrorFile = 'errors_bulk')
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
    monosDir = os.listdir(monolayer_path)

    if bulks in os.listdir('inputs'):
        bulksDir = os.listdir('inputs/bulks')
    else:
        bulks = []
    monos = {}
    bulks = {}
    # making direcotries of types of chem sub directories 
    os.mkdir('hull')
    os.mkdir('onRatio')
    # read in monolayer and bulk motifs
    for mon in monosDir:
        monos[mon] = Structure.from_file('inputs/monolayers/'+mon)
    for bulk in bulksDir:
        bulks[bulk] = Structure.from_file('inputs/bulks/'+mon)
    # perform the chem sub based on the input X and M lists 
    for monolayer in monosDir:
        os.mkdir(monolayer)
        os.chdir(monolayer)
        for M in M_elements:
            for X in X_elements:
                coreElements = [M,X]
                if ratio[0] == 1 and ratio[1] == 1:
                    os.mkdir(M+X)
                    os.chdir(M+X)
                elif ratio[0] == 1:
                    os.mkdir(M+X+str(ratio[1]))
                    os.chdir(M+X+str(ratio[1]))
                elif ratio[1] === 1:
                    os.mkdir(M+str(ratio[0])+X)
                    os.chdir(M+str(ratio[0])+X)
                else:
                    os.mkdir(M+str(ratio[0])+X+str(ratio[1]))
                    os.chdir(M+str(ratio[0])+X+str(ratio[1]))
                struct = monos[monolayer]
                spec = struct.species
                unique = []
                for specie in spec:
                    if specie not in unique:
                        unique.append(specie)
                if len(unique) != len(coreElements:
                    print monolayer+' DOES NOT HAVE '+len(coreElements)+' ELEMENTS. ENDING LOOP'
                    break
                for element in unique:
                    struct.replace({element:coreElements[unique.index(element)]})
                struct.to('POSCAR','POSCAR')
                write_potcar(pathToPOTCARS)
                relax(submit=False)
                competing = get_competing_phases()
                os.chdir('../../')
                ## test for on Ratio : if competing phases is just 1 
                if len(competing) == 1:
                    os.chdir('onRatio')
                    newStruct = get_struct_from_mp(competing[0][[1]) # assume the  
                                                                     # lowest hull distance 
                    os.mkdir(competing[0][0])    
                    os.chdir(competing[0][0])
                    newStruct.to('POSCAR','POSCAR')  # an output on Ratio
                                                     # chem sub 
                    # write_potcar(pathToPOTCARS) no need to run/write_potcar
                    # relax(submit=False)
                    os.chdir('../')
               ## how about more than 2 ? this would just be the case in ternaries, etc.
               elif len(competing) == 2:
                    os.chdir('hull')
                    newStruct1 = get_struct_from_mp(competing[0][[1])
                    os.mkdir(competing[0][0])
                    os.chdir(competing[0][0])
                    newStruct1.to('POSCAR','POSCAR')
                    write_potcar(pathToPOTCARS)
                    relax(submit=False)
                    os.chdir('../')
                    newStruct2 = get_struct_from_mp(competing[1][[1])
                    os.mkdir(competing[1][0])
                    os.chdir(competing[1][0])
                    newStruct2.to('POSCAR','POSCAR')
                    # write_potcar(pathToPOTCAR)
                    # relax(submit=False)
                    os.chdir('../')
                os.chdir('../')
    for bulk in bulksDir:
        os.mkdir(bulk+'_bulk')
        os.chdir(bulk)
        for M in M_elements:
            for X in X_elements:
                coreElements = [M,X]
                if ratio[0] == 1 and ratio[1] == 1:
                    os.mkdir(M+X)
                    os.chdir(M+X)
                elif ratio[0] == 1:
                    os.mkdir(M+X+str(ratio[1]))
                    os.chdir(M+X+str(ratio[1]))
                elif ratio[1] === 1:
                    os.mkdir(M+str(ratio[0])+X)
                    os.chdir(M+str(ratio[0])+X)
                else:
                    os.mkdir(M+str(ratio[0])+X+str(ratio[1]))
                    os.chdir(M+str(ratio[0])+X+str(ratio[1]))
                struct = monos[monolayer]
                spec = struct.species
                unique = []
                for specie in spec:
                    if specie not in unique:
                        unique.append(specie)
                if len(unique) != len(coreElements:
                    print monolayer+'_bulk DOES NOT HAVE '+len(coreElements)+' ELEMENTS. ENDING LOOP'
                    break
                for element in unique:
                    struct.replace({element:coreElements[unique.index(element)]})
                struct.to('POSCAR','POSCAR')
                # write_potcar(pathToPOTCARS)
                # relax(submit=False)
                os.chdir('../../')


if __name__=='__main__':
    Ms = ['Ir']
    Xs = ['I','Cl']
    ratio = [1,2]
    makePoscars(Ms, Xs, ratio, api_key)




