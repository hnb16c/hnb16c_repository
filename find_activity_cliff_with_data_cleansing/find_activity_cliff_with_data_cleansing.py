# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 16:09:16 2018

Data Cleansing 
Input:
    
(CMPD_CHEMBLID,CANONICAL_SMILES,pKi) 
CHEMBL230981,[O-][n+]1ccccc1CC[C@@H](NS(=O)(=O)Cc2ccccc2)C(=O)NCC(=O)NCc3cccc(Cl)c3,5
CHEMBL340577,COc1ccc(cc1)n2nnnc2C(=O)Nc3ccc(cc3)c4ccccc4S(=O)(=O)N,5
CHEMBL62509,C[C@@H]1CC[C@@H](CC1)Oc2cccc(c2O)c3nc4cc(C(=N)N)c(F)cc4[nH]3,5
CHEMBL316311,CCSCC1(CC(=NO1)c2cccc(c2)C(=N)N)C(=O)Nc3ncc(cn3)c4ccccc4S(=O)(=O)N,5


@author: Hiroshi Nakano
"""

import sys
import numpy as np
import rfrag2,indexing2,post01_make_exclusive_MMP02,post02_add_pKi_to_MMP02

def heavy_atom_count(smi):

  m = Chem.MolFromSmiles(smi)
  return m.GetNumAtoms()


def mergemultidata(file0,file1):

    smidict={}
    pKidict={}
    f0=open(file0)
    line=f0.readline()
    #print(line)
    for line in f0:
        line=line.rstrip()
        #print(line)
        id, smi, pKi = line.split(',')
        smidict.setdefault(smi,[]).append(id)
        pKidict[id]=pKi
    f0.close()

    f1=open(file1,'w')
    for uniqsmi in smidict:
        idvec = smidict[uniqsmi]
        if len(idvec) == 1:
            line  = idvec[0] + ','
            line += uniqsmi +','
            line += pKidict[idvec[0]]+'\n'
            f1.write(line)
        elif len(idvec) > 0:
            nmulti = len(idvec)
            pKivec = np.zeros(nmulti,'float')
            #print(idvec)
            for imulti in np.arange(nmulti):
                #print([idvec[imulti]])
                pKivec[imulti]=float(pKidict[idvec[imulti]])
            #print(pKivec)
            if np.max(pKivec) - np.min(pKivec) <= 1:
                unifiedpKi=np.sqrt(np.mean(pKivec**2))
                #print(unifiedpKi)
            #else:
                #print('DISCARDED')
            #print(smidict[uniqsmi])
            line  = idvec[0] + ','
            line += uniqsmi +','
            line += str(unifiedpKi)+'\n'
            f1.write(line)
    f1.close()
    
if __name__ == '__main__':
    #file0 = 'a10-thrombin-ChEMBL_24_bioactivity-18_11_01_08.csv'
    file0 = 'sample_dataset.csv'
    file1 = 'st01'+file0
    file2 = 'st02'+file0
    file3 = 'st03'+file0
    file4 = 'st04'+file0
    file5 = 'st05'+file0
    
    ## Merge the multple data with the same SMILES.
    print("merge multi data")
    mergemultidata(file0,file1)

    ## fragmentation with rfrag2.py 
    print("rfrag2")
    rfrag2.main(file1,file2)
    
    ## make MMP
    print("indexing2")
    indexing2.main(file2,file3)
    
    ## merge MMPs consit of the same molecules-pairs into the largest shared fragment
    print("make exclusive MMP02")
    post01_make_exclusive_MMP02.main(file3,file4)
    
    ## add pKi values to the exclusive MMPs
    print("add pKi to MMP02")
    post02_add_pKi_to_MMP02.main(file1,file4,file5)
    
    print("Data cleansing finished.")
    
    
