## Beta chemical substitution module 
## that can be used to perform chemical substitution 
## using pymatgen tools and names the file according to 
## the comment line provided in a structural motif 


from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from glob import glob
import os

def chem_sub_motif(motif_files, substituents, sub_both=False):
   """
   chemical substitution of a list of motifs 
   with options for the subsituents
   Args:
       motif_files (list): list of poscar file paths 
       substituents (dict): Eg: for a motif AB this dict 
                           should be 
                           {'Mo':['Ta', 'V'], 
                            'N': ['N', 'Cl'] }
       NOTE: all motifs should be of same number of components,
             and should have the same species 
             that is works if all motifs are binaries, ternaries, etc.
             of the same components 
             not a mix of these 
       sub_both (bool) : TODO if necessary, applicable for cases like 
                         CdCl2 --> MnF2 in one go 
   Returns:
       List of Poscar objects  
   """
   print 'substituents {0}'.format(substituents)
   chem_sub_poscars = [] 
   for motif in motif_files:
#      if sub_both:
#          p = Poscar.from_file(motif)
#          subp = substituents.keys()
#          p.structure.replace_species({Element(subp[0]):Element(subs_posit)})
      for component in substituents.keys():
          for s in substituents[component]:
              p1 = Poscar.from_file(motif)
              print 'substitute {0} for {1}'.format(component, s)
              p1.structure.replace_species({Element(component):Element(s)})
              p1.comment = '_'.join([str(p1.structure.formula).replace(' ', '_'), p1.comment])
              print 'created {}'.format(p1.comment)
              chem_sub_poscars.append(p1)

   print len(chem_sub_poscars)
   return chem_sub_poscars

if __name__=='__main__':
   
   motif_files = ['1T.vasp', '2H.vasp']
   substituents = {'Mo':['Fe', 'Ta'], 'N': ['Cl', 'F'] }
   poscars = chem_sub_motif(motif_files, substituents)

   #print poscars


#for f in glob.glob('*.vasp'):
#   for v in ['Ta', 'Y']:
#      p = Poscar.from_file(f)
#      p.structure.replace_species({Element('Mo'):Element(v)})  
#      p.comment = v+'N2'+'_'+p.comment
#      p.write_file(v+'N2'+os.sep+p.comment+".vasp")

