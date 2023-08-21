# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:01:54 2022
@author: ZR48SA
"""

#%% change directory to script directory
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()))
print(os.getcwd())
from Download import *
from setup_NCBI_taxonomy import *
from prep_Database import *
from make_DIAMOND_database import *

#%% Step 1: Download DIAMOND

#diamond_path=download_diamond(path=False)
diamond_path=r"C:\Users\LocalAdmin\Desktop\Spyder\NovoLign\prefinal_NovoLign_code_v19_HK02_MP\Setup\diamond\diamond\\"

#%% Step 2: Setup NCBI taxonomy

# names,nodes=download_ncbi_taxdump(path=False)
# ncbi_taxonomy_path=parse_NCBI_taxonomy(names,nodes)
ncbi_taxonomy_path=r"C:\Users\LocalAdmin\Desktop\Spyder\NovoLign\prefinal_NovoLign_code_v19_HK02_MP\Setup\ncbi_taxonomy\parsed_ncbi_taxonomy.tsv"
    
#%% Step 3: Setup Diamond database

#Database choices: "RefSeq", "NCBI_NR", "Swiss-Prot", "TrEMBL", "UniProt", "UniRef100", "UniRef90", "UniRef50"
#fasta_files=download_db("Swiss-Prot")

#merged_database=merge_files(fasta_files,delete_old=False)
#merged_database=r"C:\Users\LocalAdmin\Desktop\Spyder\NovoLign\prefinal_NovoLign_code_v19_HK02_MP\Setup\Swiss-Prot\uniprot_sprot.fasta"
#prepped_database=prep_database(merged_database,delete_old=False)
prepped_database=r"C:\Users\LocalAdmin\Desktop\Spyder\NovoLign\prefinal_NovoLign_code_v19_HK02_MP\Setup\SwissProt\uniprot_sprot_NoAmb_IJeqL.fasta"
outpath=r"C:\Users\LocalAdmin\Desktop\Spyder\NovoLign\prefinal_NovoLign_code_v19_HK02_MP\Setup\SwissProt\diamond_db"
diamond_database_path=make_diamond_database(diamond_path,prepped_database,outpath,delete_old=False)

#%% Step 4: Write filepaths to file

# path_df=pd.DataFrame([["diamond_path",                  diamond_path],
#                       ["ncbi_taxonomy_path",      ncbi_taxonomy_path],
#                       ["diamond_database_path",diamond_database_path]],columns=["name","path"]).set_index("name")

# path_df.to_csv(Path(basedir,"setup_filepaths.tsv"),sep="\t")

print("done")