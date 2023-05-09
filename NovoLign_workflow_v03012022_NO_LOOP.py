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



    
# =============================================================================
# BEGIN DEFINE PARAMETERS
# =============================================================================
    

##### Use default setup #####
Default=True           #if True, will look for setup file and overwrite manual filepaths.  
Write_to_database=True #if True, will add and write 

##### Input manual filepaths ####
#(Set Default to False)
diamond_path=         str(Path(basedir,"Setup","diamond","diamond"))                             # placeholder path to diamond executable
ncbi_taxonomy_path=   str(Path(basedir,"Setup","ncbi_taxonomy","parsed_ncbi_taxonomy.tsv"))      # placeholder path to parsed ncbi taxonomy
diamond_database_path=str(Path(basedir,"Setup","Swiss-Prot","Swiss-Prot.dmnd"))                  # placeholder path to diamond database
    
 
##### Performance parameters ####
PEAKS_Score_cuttoff=50    # minimum ALC(%)
bit=30                    # minimum bitscore


Temporary_directory=basedir #location where DIAMOND writes temmporary indices to
input_files=glob.glob("".join((basedir,"\Input_*"))) # location of input folders

# =============================================================================
# END DEFINE PARAMETERS
# =============================================================================    
    
#load NCBI parsed taxonomy
ranks=["superkingdom","phylum","class","order","family","genus","species"]
ncbi_taxdf=pd.read_csv(ncbi_taxonomy_path,sep="\t").fillna("")

#default diamond output columns
output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]
if Write_to_database: output_columns+=["full_sseq"]

#attempt to read filepaths from Setup
if Default:
    filepaths_df_path=str(Path(basedir,"setup_filepaths.tsv"))
    if os.path.exists(filepaths_df_path):
        filepaths_df=pd.read_csv(filepaths_df_path,sep="\t").set_index("name")
        diamond_path         =filepaths_df.loc["diamond_path","path"]
        ncbi_taxonomy_path   =filepaths_df.loc["ncbi_taxonomy_path","path"]
        diamond_database_path=filepaths_df.loc["diamond_database_path","path"]
        
    else:
       print("Setup filepaths not found, please execute run_Setup.py in the folder Setup, with your database of choice")
       print("Or copy-paste absolute paths to diamond executable, diamond database, and parsed ncbi taxonomy file manually in the 3 lines below:")


# =============================================================================
# START NOVOLIGN PIPELINE
# =============================================================================

from write_to_fasta import write_to_fasta
from diamond_alignment import align
from process_alignment import Process_diamond_alignment, lca
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
    fasta_file,denovo_peptides=write_to_fasta(de_novo_file,database_searching_file,Output_directory,PEAKS_Score_cuttoff=PEAKS_Score_cuttoff,add_decoy=True)  #Rewrite with key arguments not positional
    alignments=align(fasta_file,Output_directory,diamond_path=diamond_path,database_path=diamond_database_path,output_columns=output_columns)
    target_decoys=Process_diamond_alignment(alignments,Output_directory,minimum_bitscore=bit,output_columns=output_columns)

    """
    2/5 LCAs for alignments
    """
    print("Step 2 of 5: Construct LCAs")
    for file in target_decoys: 
        # 1. Conventional lca
        denovo_peptides_lca=lca(file,Output_directory,denovo_peptides=denovo_peptides,method="standard")
        # 2. Bitsore lca 
        denovo_peptides_blca=lca(file,Output_directory,denovo_peptides=denovo_peptides,method="focused",weight_column="bitscore")
        # # 3. Weighted lca
        denovo_peptides_wlca=lca(file,Output_directory,denovo_peptides=denovo_peptides,method="weighted",weight_column="weights") 
        # # 4. Decoy corrected weighted lca
        denovo_peptides_cwlca=lca(file,Output_directory,denovo_peptides=denovo_peptides,method="weighted",weight_column="corrected_weights") 

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
