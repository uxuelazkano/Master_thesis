#***********************************************
# TO USE THIS SCRIPT:
# module load RDKit/2017.03.2-Python-3.5.2
# module load Python/3.5.2
#***********************************************


import pandas as pd
import numpy as np
import gzip, os
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold # scaffold calculator
from rdkit.Chem import Recap # fragmentator
from rdkit.Chem import PandasTools # molecules in df viewer
from rdkit.Chem import Draw # molecules drawer
from rdkit.Chem import rdFMCS # maximum common substructure found
from rdkit.Chem import rdMolDescriptors # count rings
from rdkit import DataStructs # similarity measure
from rdkit.Chem.Fingerprints import FingerprintMols # similarity measure
from rdkit.Chem import Descriptors # to calculate MW
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm 
import itertools
import csv
import sys
import datetime
import pickle
import subprocess
from cycler import cycler
import pymysql
from itertools import *
import ast
import seaborn as sns

class MyUniprot:
    """ 
    Preparation of input file (Uniprot, Smiles) for fragmentation for each Uniprot accession.
    """

    def __init__(self, uniprot):
        self.uniprot = uniprot

    
    def get_smiles(self):
        SELECT = ["compound_structures.canonical_smiles"]
        FROM = ["component_sequences", "target_components","target_dictionary","assays","activities","molecule_dictionary","compound_structures","compound_records", "molecule_hierarchy"]
        WHERE = ["(component_sequences.accession LIKE '{p1}%')".format(p1=self.uniprot),
        "component_sequences.component_id=target_components.component_id","target_components.tid=target_dictionary.tid","target_dictionary.tid=assays.tid",
        "assays.assay_id=activities.assay_id","activities.record_id=compound_records.record_id","compound_records.molregno=molecule_dictionary.molregno",
        "compound_records.molregno=molecule_hierarchy.molregno", "molecule_hierarchy.parent_molregno=compound_structures.molregno","activities.pchembl_value IS NOT NULL",
        "(activities.standard_type='IC50' OR activities.standard_type='EC50' OR activities.standard_type='Ki' OR activities.standard_type='Kd')"]
        
        sql = "SELECT DISTINCT {s} FROM {f} WHERE {w};".format(s=", ".join(SELECT), # repeat query for patent_id w/o hypen
                                                      f=", ".join(FROM),
                                                      w=" AND ".join(WHERE))
        print(sql)
        if cursor.execute(sql) != 0:
            smiles = [pt["canonical_smiles"] for pt in cursor]
            return(smiles)

        else: # empty set
            return([])




class MySmiles:
    """ 
    SMILES for fragmentation.
    """

    def __init__(self, smiles):
        self.smiles = smiles

    def get_mol(self):
        mol = Chem.MolFromSmiles(self.smiles)
        return(mol)
    
    def is_mol(self): # 1st requirement
        if self.get_mol() is None: # if syntax errors
            return(False)
        else:
            return(True)
            
    def is_small(self):
        if self.is_mol() is False: # 2nd requirement
            return(False)
        elif Descriptors.MolWt(self.get_mol()) <= float(600):
            return(True)
        else:
            return(False)

    def has_2_rings(self): # 3th requirement
        if self.is_mol() is False:
            return(False)
        elif Chem.GetSSSR(self.get_mol()) < 2: # drug have normally a min. of 2 rings
            return(False)
        else:
            return(True)
        
    # -- requirements for fragmentation step finished
    
    def get_fragments(self):
        fragments = None
        if False not in [self.is_mol(), self.is_small(), self.has_2_rings()]: # 3 requirements fulfilled
            fragments = []          
            
            # 1st add scf of the fragments
            hierarch = Recap.RecapDecompose(self.get_mol())
            ks = hierarch.children
            for s, obj in ks.items():
                m = obj.mol
                if (m is None) or (Chem.GetSSSR(m) < 2): 
                    continue
                # Fragments into scaffolds conversion
                try:
                    core = MurckoScaffold.GetScaffoldForMol(m)
                except ValueError: # scf calculation not possible
                    continue
                smiles_scf = Chem.MolToSmiles(core)
                if Chem.MolFromSmiles(smiles_scf) is None: 
                    continue
                fragments.append(smiles_scf)
                
            # 2nd add scf of itself
            try:
                core = MurckoScaffold.GetScaffoldForMol(self.get_mol())
                smiles_scf = Chem.MolToSmiles(core)
                if Chem.MolFromSmiles(smiles_scf) is not None:
                    fragments.append(smiles_scf)
            except ValueError: # scf calculation not possible
                pass                
            
        return(fragments)    
    
    def get_murckoscf(self):
        core = None
        if False not in [self.is_mol(), self.is_small(), self.has_2_rings()]: # 3 requirements fulfilled
            try:
                core = MurckoScaffold.GetScaffoldForMol(self.get_mol())
            except ValueError: # scf calculation not possible
                core = None
            if core is not None:
                core = Chem.MolToSmiles(core)       
        return(core)

    def has_substructure(self, mol_substructure):
        match = None
        if None not in [self.get_mol(), mol_substructure]:
            match = self.get_mol().HasSubstructMatch(mol_substructure)
        return(match)
    
    def get_canonical_smiles(self):
        canonical_smiles = None
        if self.is_mol() is not False:
            canonical_smiles = Chem.MolToSmiles(self.get_mol())
        return(canonical_smiles)         


# Open Chembl mysql connection for checking part
connection = pymysql.connect(host='blackadder', database='chembl_23', user='syspharm')
cursor = connection.cursor(pymysql.cursors.DictCursor)

###FRAGMENTATION
# Run scaffold extraction structuresprotocol for each Uniprot (in for-loop)
u = MyUniprot('P03372')
smiles = u.get_smiles()
#print(smiles)
print(u, len(smiles))


dfSmiles = pd.DataFrame({}, columns=["UNIPROT", "MOL", "FRAGMENTS"])
dfOut = pd.DataFrame({}, columns=["UNIPROT", "FRAGMENTS", "FREQ"])
print("Get Fragments")
fragments_dic = {} #Para hacer el count y crear el threshold
FragmnetSmiles = {}
#Convertir smiles en canonical de RDkit y hacer un set 
smiles = set( [ MySmiles(s).get_canonical_smiles() for s in smiles if MySmiles(s).is_mol() ] )
print(len(smiles))

counter=0
counter1=0

for s in smiles:
    counter+=1
    #print(" first smiles  " + s)
    s = MySmiles(s)
    fragments = s.get_fragments()
    if (fragments is not None) and (len(fragments)>0) :
        for sc in fragments:
            if sc not in FragmnetSmiles.keys():
                counter1+=1
                print
                FragmnetSmiles[sc]=[]
            FragmnetSmiles[sc].append(s.smiles)
            if sc not in fragments_dic.keys():
                fragments_dic[sc] = 0
            fragments_dic[sc] += 1
    continue

#thr = np.amax(list(fragments_dic.values()))/3
#high_fragments = [ f for f, freq in fragments_dic.items() if freq >= thr ]
#print(high_fragments)  
print("Counter is, smile number: " + str(counter))
print("Counter1 is, number of fragments in smile: "+ str(counter1))
print("Number of molecules in the analysis of Q12809 is:  " + str(len(smiles)))

for key, values in FragmnetSmiles.items():
    for val in values:
        dfSmiles = dfSmiles.append({'UNIPROT':u.uniprot,'MOL':val, 'FRAGMENTS':key}, ignore_index=True)



for key, value in fragments_dic.items():
    dfOut = dfOut.append({'UNIPROT':u.uniprot,'FRAGMENTS':key, "FREQ":value}, ignore_index=True)

#for f in high_fragments: # to visualize 2D structures
#    df_2D = df_2D.append({"name": "high_fragment_sure", "smiles": f}, ignore_index=True)
#scf_plot(p, fragments_dic)



###CLUSTERING SCAFFOLDS
#Find substructures and group them in a dictionary

FragmentsLista=set(dfOut.FRAGMENTS)
FragmentsMol=[Chem.MolFromSmiles(s) for s in FragmentsLista]
match_dictionary={}

print("Start substructure finding: " + str(datetime.datetime.now()))
for m, n in permutations(FragmentsMol, 2):
    mm=Chem.MolToSmiles(m)
    nn=Chem.MolToSmiles(n)
    if nn not in match_dictionary.keys():
        match_dictionary[nn]=[]
    if m.HasSubstructMatch(n):
        match_dictionary[nn].append(mm)

print("Finish substructure finding: "+ str(datetime.datetime.now()))
print("Start similarity: " + str(datetime.datetime.now()))

#Remove keys with empty lists
match_dictionary={k:v for k,v in match_dictionary.items() if len(v)>0}

#Remove subgroups an leave only the big group with the smallest smile as key
remove=[]
base_keys = match_dictionary.keys()
for key, values in match_dictionary.items():
    for words in values:
        for base_key in base_keys:
            if base_key in words:
                remove.append(base_key)

for i in remove:
    if i in match_dictionary.keys():
        del match_dictionary[i]

#Create new df with the smiles and corresponding scaffold for each one
df1=pd.DataFrame({},columns=['UNIPROT','FREQ','FRAGMENTS','SCAFFOLD'])
print("Creating scaffolds")
for i, row in dfOut.iterrows():
    sc=[]
    for key, values in match_dictionary.items():
        if row['FRAGMENTS'] in values:
            sc.append(key)
    if len(sc)==0:
        sc.append(row['FRAGMENTS'])
    for s in sc:
        df1=df1.append({'UNIPROT':row['UNIPROT'],'FREQ':row['FREQ'], 'FRAGMENTS':row['FRAGMENTS'], 'SCAFFOLD':s},ignore_index=True)


dfSmiles.drop(dfSmiles.columns[[0]],axis=1, inplace=True)

df1.sort_values("SCAFFOLD")

print(dfSmiles.head(10))
print(df1.head(10))

M = len(dfSmiles['MOL'].unique().tolist())
F = len(dfSmiles.FRAGMENTS.tolist())
print("Total nº of frags is: "+str(F))
print("Total nº of molecules in the target: " + str(M))

dictionary={}
for sc in set(df1.SCAFFOLD.tolist()):
    dictionary[sc]={}
    mols = [] 
    frag = df1[df1.SCAFFOLD == sc]["FRAGMENTS"].tolist() 
    for f in frag:
        mols += dfSmiles[dfSmiles.FRAGMENTS == f]["MOL"].tolist()

    #Fill dict
    dictionary[sc]={"frag": (len(mols)*100)/F, "mols": mols}
    
sorted_scaffold = sorted(dictionary, key=lambda x: dictionary[x]["frag"], reverse=True)

print("Start similarity calc: "+ str(datetime.datetime.now()))

dictionary1={}
for key, values in dictionary.items():
    if key not in dictionary1.keys():
        dictionary1[key]={}
    mols=list(dictionary[key]['mols'])
    sim=[]
    
    for m, n in combinations(mols, 2):
        mm = Chem.MolFromSmiles(m)
        nn = Chem.MolFromSmiles(n)
        fp1 = AllChem.GetMorganFingerprint(mm,2)
        fp2 = AllChem.GetMorganFingerprint(nn,2)
        similarity =DataStructs.DiceSimilarity(fp1,fp2)
        sim.append(similarity)
        dictionary1[key]=sim

print(len(dictionary1))



#dictionary1[key]=sim
print("Finish similarity calc: "+ str(datetime.datetime.now()))


sorted_scaffold = sorted(dictionary1, key=lambda k: len(dictionary1[k]), reverse=True)
print(sorted_scaffold)

counter=0

for i in sorted_scaffold:
    f, axes = plt.subplots(2, 2, figsize=(7, 7), sharex=True)
    sns.despine(left=True)
    print(i)
    counter+=1
    #plt.hist(dictionary1[i],bins='auto', histtype='stepfilled', label=i)
    d=dictionary1[i]
    #sns.distplot(d, color="m")
    #plt.ylabel ('%mol')
    #plt.xlabel ('Similarity')
    #plt.legend (bbox_to_anchor=(1, 1), loc="upper right", borderaxespad=0.)
    #plt.show()
    if counter==1:
        sns.distplot(d, hist=False, color="m", kde_kws={"shade": True}, ax=axes[0, 0])
    if counter==2:
        sns.distplot(d, hist=False, color="m", kde_kws={"shade": True}, ax=axes[0, 1])
    if counter==3:
        sns.distplot(d, hist=False, color="m", kde_kws={"shade": True}, ax=axes[1, 0])
    if counter==4:
        sns.distplot(d, hist=False, color="m", kde_kws={"shade": True}, ax=axes[1, 1])
    if counter ==5:
        break


# Close connection 
cursor.close()
connection.close()
