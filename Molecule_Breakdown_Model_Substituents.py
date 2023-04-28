#!/usr/bin/env python
# coding: utf-8

# # **_Substituent Database Building_**

# ## Environment Setting Up and Packages Importing:
#         !conda activate fragment-env 
#         Please create this environment from the .yml file: fragment.yml

# In[ ]:


import sys
import time

import pandas as pd
from tqdm import tqdm 
import pickle
import numpy as np
from collections import Counter

#RDKit:
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit import Chem
from rdkit.Chem.rdmolops import *
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import MolsToGridImage


# ### 1. Input the Data (usually for multiple-input-file cases):

# In[ ]:


# To check the time running the code:
start_time = time.time()

scr_path = sys.argv[1]
filename = (scr_path.split('/')[-1]).split('.')[0]

df_input = pd.read_csv(scr_path,  sep = ' ', header = None, names = ['SMILES'])



# Remove the duplicates in the original dataset:
df_input.drop_duplicates(subset = "SMILES", keep = 'last', inplace = True)
df_input = df_input.reset_index(drop = True)
df_input.head(3)



# To check the size of the database:
len(df_input)


# In[ ]:


output1 = 'This databse', filename, 'contains', len(df_input), 'compounds'


# ### 2. Create a Data Output Dataframe:

# In[ ]:


df_output = pd.DataFrame({'Input_Smiles':df_input['SMILES']}) 
df_output['Fragment_Distribution'] = 'fragment_distribution'
df_output


# ### 3. Break Molecules from All SMILES:
# 

# #### (1) Remove All Bonds Between a Cyclic Atom and an Acyclic Atom:
#         To obtain all the substituents on rings
# #### (2) Remove All Bonds Between Any Two Cyclic Atoms:
#         To break the ring structures
#        

# In[ ]:


# Neutralizing Molecules:
def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


# In[ ]:


count = 0
my_list_invaild = []

for q in tqdm(range(0, len(df_output)), desc = 'Loop 1'):

    try:  
        # 1.remove all the non-bond ions:
        qsmiles = df_output['Input_Smiles'][q]
        qsmiles = qsmiles.replace('.[Li+]', '')
        qsmiles = qsmiles.replace('.[Na+]', '')
        qsmiles = qsmiles.replace('.[K+]', '')
        qsmiles = qsmiles.replace('.[F-]', '')
        qsmiles = qsmiles.replace('.[Cl-]', '')
        qsmiles = qsmiles.replace('.[Br-]', '')
        qsmiles = qsmiles.replace('.[I-]', '')

        qmol = Chem.MolFromSmiles(qsmiles)
        neutralize_atoms(qmol)
        qmolw = Chem.RWMol(qmol)
 
        # 2.make sure there is no more dissociation 
        if '.' not in qsmiles:

            if qmol and qmolw:

                # 3. Save all the cyclic atom in a list
                cyclic_list=[]
                for atom in qmol.GetAtoms():
                    if atom.IsInRing() == True:
                        i=atom.GetIdx()
                        c = qmol.GetAtomWithIdx(i).GetSymbol()
                        cyclic_list.append(c)
    
                # 4.Remove all the bonds between two cyclic atoms in the same ring:     
                for bond in qmol.GetBonds():
                    
                    if bond.IsInRing() == True:
                        qmol_begin=bond.GetBeginAtomIdx()
                        qmol_end=bond.GetEndAtomIdx()
                        bondatom=(qmol_begin, qmol_end)
                        bond_a=qmol.GetAtomWithIdx(qmol_begin).GetSymbol()
                        bond_b=qmol.GetAtomWithIdx(qmol_end).GetSymbol()
                        if qmol.GetAtomWithIdx(qmol_begin).IsInRing() == qmol.GetAtomWithIdx(qmol_end).IsInRing() == True:
                            qmolw.RemoveBond(qmol_begin, qmol_end)
                                    
                    else:
                        qmol_begin=bond.GetBeginAtomIdx()
                        qmol_end=bond.GetEndAtomIdx()
                        bondatom=(qmol_begin, qmol_end)

                        # 5.Remove the cyclic atom in the subtituents:
                        if  qmol.GetAtomWithIdx(qmol_begin).IsInRing() != qmol.GetAtomWithIdx(qmol_end).IsInRing():
                            qmolw.RemoveBond(qmol_begin, qmol_end)

                        # 6.Remove the bonds between two cyclic atons in different rings:
                        elif qmol.GetAtomWithIdx(qmol_begin).IsInRing() == qmol.GetAtomWithIdx(qmol_end).IsInRing() == True:
                            qmolw.RemoveBond(qmol_begin, qmol_end)
                            
                new_mol=qmolw.GetMol()
                new_smiles = Chem.MolToSmiles(new_mol)
                frag_simles = new_smiles.split('.')
                mylist = frag_simles

                mylist=list(map(lambda item: item.replace('[C@@H]','C'), mylist))
                mylist=list(map(lambda item: item.replace('[C@H]','C'), mylist))
                mylist=list(map(lambda item: item.replace('[C@@]','C'), mylist))
                mylist=list(map(lambda item: item.replace('[C@]','C'), mylist))
                mylist=list(map(lambda item: item.replace('[P@@H]','P'), mylist))
                mylist=list(map(lambda item: item.replace('[P@@]','P'), mylist))
                mylist=list(map(lambda item: item.replace('[P@]','P'), mylist))
                mylist=list(map(lambda item: item.replace('[S@@H]','S'), mylist))
                mylist=list(map(lambda item: item.replace('[S@@]','S'), mylist))
                mylist=list(map(lambda item: item.replace('[S@]','S'), mylist))
                mylist=list(map(lambda item: item.replace('[N@]','N'), mylist))
        
                mylist=list(map(lambda item: item.replace('[nH3]','N'), mylist))

                mylist=list(map(lambda item: item.replace('C-C','CC'), mylist))
                
                mylist=list(map(lambda item: item.replace('c','C'), mylist))
                mylist=list(map(lambda item: item.replace('n','N'), mylist))
                mylist=list(map(lambda item: item.replace('s','S'), mylist))
                mylist=list(map(lambda item: item.replace('o','O'), mylist))
                mylist=list(map(lambda item: item.replace('p','P'), mylist))

                mylist=list(map(lambda item: item.replace('[ZN','[Zn'), mylist))
                mylist=list(map(lambda item: item.replace('[MO','[Mo'), mylist))
                mylist=list(map(lambda item: item.replace('[AS','[As'), mylist))
                mylist=list(map(lambda item: item.replace('[AC','[Ac'), mylist))
                mylist=list(map(lambda item: item.replace('[SN','[Sn'), mylist))
                mylist=list(map(lambda item: item.replace('[IN','[In'), mylist))
                mylist=list(map(lambda item: item.replace('[te','[Te'), mylist))

                c1 = Counter(mylist)
                c2 = Counter(cyclic_list)
                mylist_new = c1-c2

                #repeat-count:
                my_dict = mylist_new

                # 7.Count each fragment only once (single-count)
                my_dict = {x: 1 for x in my_dict}

                key_list=list(my_dict.keys())
                count_list=list(my_dict.values())
                df_output['Fragment_Distribution'][q]= my_dict

        else:
            count += 1
            my_list_invaild.append(q)
            continue

    except:
        count += 1
        my_list_invaild.append(q)
        continue

print(count, 'invalids')


# ### 4. Wrangling the Data Obtained :
# 

# In[ ]:


# 1.Remove the invalid data: 
new_df = df_output.drop(my_list_invaild)

# Reset the index:
new_df = new_df.reset_index(drop = True)

print('After remove invalid data:', len(new_df))


# In[ ]:


# 2.Save the pickle file:
with open (filename + '-each-cpd-breakdown-substituent-without-ROMol.pickle','wb') as f:
    pickle.dump(new_df,f)


# In[ ]:


# 3.Sort the data according to the fragment frequency:
fragment_distribution = {
    k:[d.get(k) for d in new_df.Fragment_Distribution if k in d]
    for k in set().union(*new_df.Fragment_Distribution)
}


# In[ ]:


df_result=pd.DataFrame({'Fragment':fragment_distribution}) 
k=list(df_result['Fragment'].keys())
arr=df_result.Fragment.values
l1=list()
for ls in arr:
    while None in ls:
        ls.remove(None)
        next
    else:
        l1.append(ls)

l2=list()
for i in l1:
    s=sum(i)
    l2.append(s)

df_result['Fragment_SMILES']=k
df_result['Frequency']=l2

df_result = df_result.sort_values(by=['Frequency'], ascending=False)
df_result = df_result.reset_index(drop=True)


# In[ ]:


df_result = df_result.drop(columns='Fragment')
df_result


# In[ ]:

output2 = 'The number of substituents obtained is:', len(df_result)

# In[ ]:


# 4.Save the pickle file:
# BE CAREFULLY EXCUTING -- WHEN SAVING AS PICKELS
with open (filename + '-substituents-sorted-without-ROMol.pickle','wb') as f:
    pickle.dump(df_result,f)


# In[ ]:


# 5.Finish all the code running and print the key outputs:
end_time = time.time()
duration = round((end_time - start_time)/60)
duration = round((end_time - start_time)/60)
output3 = 'This programme took', duration, 'minutes'

print(output1, output2, output3)


# ### 5. Display the Ring Fragments (ony for single-input-file cases, for multiple-input-file cases please use the 'merge' code):
#         
# 

# In[ ]:


# 1.Add ROMol to each ring fragment:
PandasTools.AddMoleculeColumnToFrame(df_result, smilesCol = "Fragment_SMILES")
df_result

# If the ROMols already have been added before, to display them again:
# PandasTools.RenderImagesInAllDataFrames(images=True)


# In[ ]:


# 2. Display the results (all or top 10000 frequent ring fragments) in a .html file
df_result_top = df_result.head(10000)

fmolport = open(filename + '-top10000-substituents-with-ROMol.html','w')
h = df_result.to_html()
fmolport.write(h)
fmolport.close()

