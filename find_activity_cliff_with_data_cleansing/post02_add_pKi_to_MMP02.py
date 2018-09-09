# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 22:36:37 2018

@author: Hiroshi Nakano

python pot02-add-pKi-to-MMP01.py < a07-mmp-results-reduced.csv > a08-mmp-results-with-pKi.csv
"""


import sys
#import re
from rdkit import Chem
#from optparse import OptionParser


def heavy_atom_count(smi):

  m = Chem.MolFromSmiles(smi)
  return m.GetNumAtoms()

#if __name__ == '__main__':
def main(filein1,filein2,fileout):
    f1=open(filein1)#'single-molucles-data.csv'

    ChEMBL_dict = {}
    iTh=0
    for iline,line in enumerate(f1):
        if iline >= 0: #Non-skip the first line
            line = line.rstrip()
            #print(iline)
            ChEMBL_id, colB, pKi = line.split(',')
            #print(ChEMBL_id)
            ChEMBL_dict[ChEMBL_id]=pKi
            if float(pKi)>=5:
                iTh+=1
    f1.close()
    f2=open(filein2)
    fout=open(fileout,'w')
    fout.write('SMILES_a,SMILES_b,id_a,id_b,transformation,shared_fragment,pKi_a,pKi_b,del_pKi,flagTh,flagAC,flagnonAC\n')

    iAC=0
    inonAC=0
    for iline,line in enumerate(f2):
        line = line.rstrip()
        core_a, core_b, id_a, id_b, smirks, context = line.split(',')
        #print(id_a)
        pKi_a = ChEMBL_dict[id_a]
        #print(id_b)
        pKi_b = ChEMBL_dict[id_b]
        del_pKi = float(pKi_a)-float(pKi_b)

        if float(pKi_a) >= 5 and float(pKi_b) >= 5: #Ki < 10uM
            flagTh = '1'
            flagAC  = (abs(del_pKi)>=2)# or (del_pKi<=-2)
            flagnonAC = (abs(del_pKi)<=1)
            if flagAC:
                flagAC='1'
                iAC+=1
            else:
                flagAC='0'
            if flagnonAC:
                flagnonAC='1'
                inonAC+=1
            else:
                flagnonAC='0'
        else:
            flagTh='0'
            flagAC='0'
        fout.write(line)
        fout.write(','+pKi_a+','+pKi_b+','+str(del_pKi)+','+flagTh+','+flagAC+','+flagnonAC+'\n')

    f2.close()
    fout.close()
    
    print('Num of MMPs with pKi larger than 5 is ' + str(iTh))
    print('Activity cliff with more than 2 digit Ki difference is ' + str(iAC))
    print('non-Activity cliff with less than 1 digit Ki difference is ' + str(inonAC))


    