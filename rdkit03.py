
# coding: utf-8

# In[20]:


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
mol1 = Chem.MolFromSmiles('C(C(=O)O)NC(=O)C(=O)O')
mol2 = Chem.MolFromSmiles('COC(=O)C(=O)NCC(=O)O')
mol3 = Chem.MolFromSmiles('C(CC(=O)O)C(C(=O)O)S')


# In[21]:


mol1=Chem.AddHs(mol1)
mol2=Chem.AddHs(mol2)
mol3=Chem.AddHs(mol3)


# In[22]:


mol1


# In[23]:


mol2


# In[24]:


mol3

