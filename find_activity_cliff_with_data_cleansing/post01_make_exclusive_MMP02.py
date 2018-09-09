# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 22:36:37 2018

@author: Hiroshi Nakano

The MMPs written by indexsing2.py are not unique.
Therefore, the program reduce the mmp-set with mutually exclusive MMPs

post01-make-exclusive-MMP02.py
"""


import sys
#import re
from rdkit import Chem
#from optparse import OptionParser


def heavy_atom_count(smi):

  m = Chem.MolFromSmiles(smi)
  return m.GetNumAtoms()


def main(filein,fileout):
    
    #mmpout = []
    mmpdict1 = {}
    mmpdict2 = {}

    #filein = sys.argv[1]
    fin = open(filein)
    fout = open(fileout,'w')
    #for iline,line in enumerate(fin):
    for iline,line in enumerate(fin):
        line = line.rstrip()
        core_a, core_b, id_a, id_b, smirks, context = line.split(',')
        
        id_mmp = id_a+id_b
        if id_mmp in mmpdict1:
            hvatom_context = heavy_atom_count(context)
            
            # switch '<' or '>'
            # '<': keep the largest context
            # '>': keep the smallest context
            if mmpdict1[id_mmp] < hvatom_context:   
                mmpdict1[id_mmp] = hvatom_context
                mmpdict2[id_mmp] = line
        else:
            mmpdict1[id_mmp] = heavy_atom_count(context)
            mmpdict2[id_mmp] = line
    #fin.close()
    #fin = open(filein)
    for line in mmpdict2.values():
        fout.write(line+'\n')
 #       if True:    
 #           print("%s,%s,%s,%s,%s,%s" %
 #                   (core_a, core_b, id_a, id_b, smirks, context))
    fin.close()
    fout.close()
