# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:55:59 2021
@author: hbckleikamp
modified: mp May 2022
further modified lcas HK 25/06/2022
created modules mp 26/06/2022
further updated code/modules mp/hk nov/dec 2022
"""

# =============================================================================
# IMPORT MODULES, DEFINE PATHS AND PARAMETERS
# =============================================================================

import glob, os
from pathlib import Path
from inspect import getsourcefile
import pandas as pd
import numpy as np
from itertools import chain
from collections import Counter
import time
from bar_graphs import graphs
start=time.time()
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
basedir=os.getcwd()


#attempt to read filepaths from Setup
filepaths_df_path=str(Path(basedir,"Setup","setup_filepaths.tsv"))
if os.path.exists(filepaths_df_path):
    filepaths_df=pd.read_csv(filepaths_df_path,sep="\t")
    diamond_path         =filepaths_df.loc["diamond_path","path"]
    ncbi_taxonomy_path   =filepaths_df.loc["ncbi_taxonomy_path","path"]
    diamond_database_path=filepaths_df.loc["diamond_database_path","path"]
    
 else:
    print("Setup filepaths not found, please execute run_Setup.py in the folder Setup, with your database of choice")
    print("Or copy-paste absolute paths to diamond executable, diamond database, and parsed ncbi taxonomy file manually in the 3 lines below:")
    
# =============================================================================
# BEGIN DEFINE PARAMETERS
# =============================================================================
    
    
    ##### Input manual filepaths ####
    diamond_path=         str(Path(basedir,"Setup","diamond"))                       #placeholder path to diamond executable
    ncbi_taxonomy_path=   str(Path(basedir,"Setup","parsed_ncbi_taxonomy.tsv"))      #placeholder path to parsed ncbi taxonomy
    diamond_database_path=str(Path(basedir,"Setup","Swiss-Prot","Swiss-Prot.dmnd"))  #placeholder path to diamond database
    
 
##### Performance parameters ####
PEAKS_Score_cuttoff=50    # minimum ALC(%)
bit=30                    # minimum bitscore
filter_dynamic_score=False
filter_alc_bit=True

Temporary_directory=basedir #location where DIAMOND writes temmporary indices to
input_files=glob.glob("".join((basedir,"\Input_*"))) # location of input folders

# =============================================================================
# END DEFINE PARAMETERS
# =============================================================================    
    
#load NCBI parsed taxonomy
ranks=["superkingdom","phylum","class","order","family","genus","species"]
ncbi_taxdf=pd.read_csv(ncbi_taxonomy_path,sep="\t")
ncbi_taxdf.columns=["OX"]+ranks+["OS"]

# =============================================================================
# START NOVOLIGN PIPELINE
# =============================================================================

from write_to_fasta import write_to_fasta
from diamond_alignment import align
from process_alignment import Process_diamond_alignment, score_lca, lca
from experiment_qc import Plot_high_scoring
# from database_qc import * | not used for parameter sweep

for infolder in input_files:
    print("analysing: "+infolder)
    
    # define input and NovoLign modules
    de_novo_file=''.join(glob.glob("".join((infolder,"\*de novo peptides.csv"))))
    database_searching_file=''.join(glob.glob("".join((infolder,"\*psm.csv"))))
    fasta_database=''.join(glob.glob("".join((infolder,"\*.fasta"))))

    """
    Prepare folders and paths
    """
    Output_directory=str(Path(Path(infolder).parents[0],Path(infolder).name.replace("Input_","Output_")))
    if not os.path.exists(Output_directory): os.mkdir(Output_directory)
    if not os.path.exists(Temporary_directory): os.mkdir(Temporary_directory)

    """
    1/5 DMD alignment
    """
    print("Step 1 of 5: Start diamond alignment")
    fasta_files,denovo_peptides=write_to_fasta(de_novo_file,database_searching_file,Output_directory,PEAKS_Score_cuttoff,add_decoy=True)
    alignments=align(fasta_files,Output_directory,DB=diamond_database_path)
    target_decoys=Process_diamond_alignment(alignments,Output_directory,bit,filter_dynamic_score,filter_alc_bit)

    """
    2/5 LCAs for alignments
    """
    print("Step 2 of 5: Construct LCAs")
    for file in target_decoys[:-1]: # ignore DN only decoy
        # 1. Conventional lca
        denovo_peptides_lca=lca(file,denovo_peptides,Output_directory)
        # 2. Bitsore lca 
        denovo_peptides_blca=score_lca(file,denovo_peptides,Output_directory,score_column="bitscore")
        # # 3. Weighted lca
        denovo_peptides_wlca=lca(file,denovo_peptides,Output_directory,weighted=True,weight_column="weights") 
        # # 4. Decoy corrected weighted lca
        denovo_peptides_cwlca=lca(file,denovo_peptides,Output_directory,weighted=True,weight_column="corrected_weights") 

    """
    3/5 Grouped taxonomy report
    """
    print("Step 3 of 5: Get DN compositions")
    files=glob.glob("".join((Output_directory,"\lca\*_combined.tsv")))
    for file in files: graphs(file,Output_directory,ranks,chain,Counter,os,np)

    """
    4/5 Check quality of fragmentation spectra
    """  
    print("Step 4 of 5: Check spectral quality")
    Plot_high_scoring(denovo_peptides,target_decoys[0],database_searching_file,Output_directory,de_novo_file)
    

    # """
    # 5/5 Check coverage obtained by database searching
    # """
    # print("Step 5 of 5: Check DB searching output")
    # Compare_Bar(denovo_peptides_blca,database_searching_file,fasta_database)

# =============================================================================
# DONE
# =============================================================================

end=time.time()
print("NovoLign completed: Run time",round((end-start)/60,2),"minutes")
globals().clear()
