#!/usr/bin/env python
# coding: utf-8

# # **_Ring Fragment Database Building_**

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
df_input


# To check the size of the database:
len(df_input)


output1 = 'This databse', filename, 'contains', len(df_input), 'compounds'



# ### 2. Create a Data Output Dataframe:

# In[ ]:


df_output = pd.DataFrame({'Input_Smiles':df_input['SMILES']}) 
df_output['Fragment_Distribution'] = 'fragment_distribution'
df_output


# ### 3. Break Molecules from All SMILES:
# 

# 
# #### (1) Remove All Bonds Between Any Two Acyclic Atoms:
#         Only keep the substituents on rings
# #### (2) Break the Ring Systems that Share One Atom:
#         Have to add the shared atom back to the sub-fragments
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

            # 3.Remove all the bonds between two acyclic atoms:
                for bond in qmol.GetBonds():
                    if bond.IsInRing() == False:
                        qmol_begin=bond.GetBeginAtomIdx()
                        qmol_end=bond.GetEndAtomIdx()
                        bondatom=(qmol_begin, qmol_end)
                        bond_a=qmol.GetAtomWithIdx(qmol_begin).GetSymbol()
                        bond_b=qmol.GetAtomWithIdx(qmol_end).GetSymbol()
                        if  qmol.GetAtomWithIdx(qmol_begin).IsInRing() == qmol.GetAtomWithIdx(qmol_end).IsInRing() == False:
                            qmolw.RemoveBond(qmol_begin, qmol_end)

                new_mol=qmolw.GetMol()
                new_smiles = Chem.MolToSmiles(new_mol)
                frag_smiles_list = []
                frag_smiles_list = frag_smiles_list + list(new_smiles.split('.'))

            # 4.Break ring systems that sharing one atom:
            result_list = []
            working_list = frag_smiles_list
            working_list_breakable_frags = []
            remove_atom_list = []
            initial_len_working_list = len(working_list)
            new_len_working_list = 0

            while new_len_working_list != initial_len_working_list:

                for f in range(1000):

                    try:
                        frag_smiles = working_list[f]
                        initial_len_working_list = len(working_list)
                            
                        if frag_smiles not in working_list_breakable_frags:
                            fragmol = Chem.MolFromSmiles(frag_smiles)
                            fragmolw = Chem.RWMol(fragmol)

                            for atom in fragmol.GetAtoms():
                                i = atom.GetIdx()
                                neighbor_list_atom = []
                                neighbor_list_idx = []
                                degree = atom.GetDegree()
                                
                                if fragmol.GetAtomWithIdx(i).IsInRing() == False and degree > 1:
                                    neighbor_list_RDKMol = fragmol.GetAtomWithIdx(i).GetNeighbors()

                                    for n in neighbor_list_RDKMol:
                                        neighbor_atom = n.GetSymbol()
                                        neighbor_idx = n.GetIdx()
                                        neighbor_list_atom.append(neighbor_atom)
                                        neighbor_list_idx.append(neighbor_idx)

                                    if all(fragmol.GetAtomWithIdx(neighbor_idx).IsInRing() == True for neighbor_idx in neighbor_list_idx) == True:

                                        working_list_breakable_frags.append(frag_smiles)

                                        for neighbor_idx in neighbor_list_idx:
                                            bond = fragmol.GetBondBetweenAtoms(neighbor_idx, i)
                                            bondtype = bond.GetBondType()
                                            bridge_atom_symbol = fragmol.GetAtomWithIdx(i).GetSymbol()
                                            fragmolw.RemoveBond(neighbor_idx, i)
                                            new_atom_idx = fragmolw.AddAtom(atom)
                                            fragmolw.AddBond(new_atom_idx, neighbor_idx, bondtype)


                                        fragmolw.RemoveAtom(i)

                                        new_fragmol = fragmolw.GetMol()
                                        new_fragsmiles = Chem.MolToSmiles(new_fragmol)

                                        working_list_split_frags = []
                                        working_list_split_frags = list(new_fragsmiles.split('.'))
     
                                        working_list = working_list + working_list_split_frags
                                        new_len_working_list = len(working_list)
                                        break

                                    else:
                                        continue
                                else:
                                    continue  
                            else:
                                new_len_working_list = len(working_list)
                        else:
                            break
                            
                    except IndexError:
                        pass
        
                c1 = Counter(working_list)
                c2 = Counter(working_list_breakable_frags)
                c3 = Counter(remove_atom_list)
                diff = c1 - c2 - c3

                result_list = list(diff.elements())
                mylist = result_list


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
                
                #repeat-count:
                my_dict = {i:mylist.count(i) for i in mylist}

                # 5.Count each fragment only once (single-count)
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
with open (filename +' -each-cpd-breakdown-ring-without-ROMol.pickle','wb') as f:
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


# 4.Divide into Ring and Non-Ring Two groups:
if_ring_check=[]  
for i in df_result['Fragment_SMILES']:
    try:
        mol = AllChem.MolFromSmiles(i)
        atoms = mol.GetAtoms()
        n_ringatoms = 0
        for atom in atoms:
            if atom.IsInRing():
                n_ringatoms += 1

            if n_ringatoms == 0:
                ringcheck = False
            else:
                ringcheck = True
        if_ring_check.append(ringcheck)
    except:
        ringcheck = False
        if_ring_check.append(ringcheck)

df_result['Ring']=if_ring_check
df_result


# In[ ]:


# 5.Obtain the results only for the ring-fragment group:
df_ring = df_result[df_result.Ring]

# Reset the index:
df_ring = df_ring.reset_index(drop=True)


# In[ ]:


output2 = 'The number of ring fragments obtained is:', len(df_ring)


# In[ ]:


# 6.Save the pickle file:
with open (filename + '-ring-fragments-sorted-without-ROMol.pickle','wb') as f:
    pickle.dump(df_ring,f)


# In[ ]:


# 7.Finish all the code running and print the key outputs:
end_time = time.time()
duration = round((end_time - start_time)/60)
output3 = 'This programme took', duration, 'minutes'

print(output1, output2, output3)



# ### 5. Display the Ring Fragments (ony for single-input-file cases, for multiple-input-file cases please use the 'merge' code):
# 

# In[ ]:


# 1.Add ROMol to each ring fragment:
PandasTools.AddMoleculeColumnToFrame(df_ring, smilesCol = "Fragment_SMILES")
df_ring

# If the ROMols already have been added before, to display them again:
# PandasTools.RenderImagesInAllDataFrames(images=True)


# In[ ]:


# 2. Display the results (all or top 10000 frequent ring fragments) in a .html file
df_ring_top = df_ring.head(10000)

fmolport = open(filename +'-top10000-ring-fragments-with-ROMol.html','w')
h = df_ring.to_html()
fmolport.write(h)
fmolport.close()

