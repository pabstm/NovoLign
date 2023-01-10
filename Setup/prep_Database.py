# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:20:18 2023

@author: hugokleikamp
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()).parents[0]) #change base directory to HybridCycler
os.chdir(basedir)
print(os.getcwd())

#%% import 


import shutil
import Bio
from Bio import SeqIO
import pandas as pd
import itertools
import random


#%%
def merge_files(folder,delete_old=False):
    

    files=[str(Path(Path(folder).parents[0],i)) for i in os.listdir(folder) if i.endswith(".fa") or i.endswith(".fasta")] 
    
    if len(files)>1:
    
        outpath=Path(folder).parents[0].name+".fa"
        
        with open(outpath,'wb') as wfd:
            for f in files:
                print(f)
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
                    
                if delete_old:
                    os.remove(f)
        
        return outpath
    
    else:
        return files[0]


def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

def prep_database(
                  Path_to_db,
                  Equate_IL=True,                                 # Equates I and J to L
                  Remove_ambiguous=["J","B","X","Z","[","(",">"], # "O","U", does not remove J in case of Equate_IL
                  Bacterial_only=False,                           # to turn on: supply here path to parsed ncbi file
                  No_Fragments=False,                             # Remove fragmented sequences  
                  Add_decoy=False,                                # append decoy of reversed or randomized peptides
                  decoy_delimiter="decoy_",
                  decoy_method="reverse",
                  Output_path=False,                              # full path to database output
                  delete_old=False,
                  
                  ):


    #for Bacterial_only (removes all viral and eukaryotic sequences from database)
    Path_to_taxonomy=Bacterial_only # path to parsed NCBI taxonomy (as created by https://github.com/hbckleikamp/NCBI2Lineage)
    if Bacterial_only:
        ranks=["superkingdom","phylum","class","order","family","genus","species"] 
        taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")
        taxdf.columns=['OX']+ranks+["OS"]
        taxa=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)["OX"].atype(int).tolist()
    
    #parse output_path
    if not Output_path:
        Output_path=Path_to_db
        if Bacterial_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
        if Remove_ambiguous: Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
        if Equate_IL:        Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
        if Add_decoy:        Output_path=Output_path.replace(".fasta","_Decoy.fasta")
    

    #read    
    recs=SeqIO.parse(Path_to_db,format="fasta")
    chunks=chunk_gen(recs)
    
    #write 
    print("writing "+Path(Output_path).stem)
    once=True
    with open(Output_path,"w") as f: #clears file if exists
        pass
    with open(Output_path,"a") as f:
    
        once=True
        for ic,c in enumerate(chunks):
            print("chunk "+ str(ic))
    
            chunk_df=pd.DataFrame([[str(r.seq),r.description] for r in c],columns=["seq","description"])

            
    
            if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
            if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Remove_ambiguous],axis=1).any(axis=1)]
            if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]

            if Bacterial_only: #only retain archaeal and bacterial sequences
                if once: #check for taxid delimiter once 
                    if chunk_df.head(100).description.str.contains("TaxID="):
                        Taxid_delimiter="TaxID="
                    else:
                        Taxid_delimiter="OX="
                    once=False
                chunk_df=chunk_df[chunk_df.description.str.contains(Taxid_delimiter)]
                chunk_df=chunk_df[chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[1]).str.split("").apply(lambda x: x[0]).astype(int).isin(taxa)]
                
                    
                
                
            
            if Add_decoy:
                decoy=chunk_df.copy()
                if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                decoy["description"]=decoy_delimiter+decoy["description"]
                chunk_df=pd.concat([chunk_df,decoy])
                
            f.write("\n"+"\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"]))

    if delete_old:
        os.remove(Path_to_db)
    return Output_path
