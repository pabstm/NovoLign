"""
module for taxonomic grouping using various LCAs
Author: mpabst
Date: 26-06-2022
"""

from __main__ import basedir, Temporary_directory
from pathlib import Path
import os
import subprocess
import shutil


def align(input_files,Output_directory, # peptides in fasta format
                      database_path,    # full path to diamond database
                      diamond_path=str(Path(basedir,"diamond")),
                      matrix_path=str(Path(basedir,"PAM70.mat")),
                      output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"],
                      select=" -k50 ", # can be --top x% or --max-target-seqs/-k
                      block_size=5,
                      index_chunks=1,
                      minimum_pident=85,
                      minimum_coverage=80,
                      minimum_bitscore=20,
                      gap_open=2,
                      gap_extend=4 
                      ):
    
    
    output_folder=str(Path(Output_directory,"diamond_alignments"))    
    if not os.path.exists(output_folder): os.mkdir(output_folder)
    
    output_files=[]
    
    if type(input_files)==str: input_files=[input_files]
    for input_file in input_files:
        output_file=str(Path(output_folder,Path(input_file).stem+".tsv"))
        
        command="cd "+'"'+basedir +'"'+ " && " + \
                "".join(['"'+diamond_path+'"',   
                " blastp -q " +'"'+input_file+'"',
                " -d "+'"'+database_path+'"',
                " -o "+'"'+output_file+'"',
                " -c" + str(index_chunks), 
                " -b "+ str(block_size),  # parameter that determines ram consumption and performance
                " "+select+" ",
                " --log ",
                " --custom-matrix "+'"'+matrix_path+'"'+" --gapopen "+str(gap_open)+" --gapextend "+str(gap_extend),  # scoring matrix path, matrix metrics
                " --algo ctg --dbsize 1 ", 
                " --id "+          str(minimum_pident),
                " --min-score "+   str(minimum_bitscore),
                " --query-cover  "+str(minimum_coverage),
                " -f 6 qseqid "+" ".join(output_columns)+" ",
                " -t "+'"'+Temporary_directory+'"'])
                
        print("diamond command:")
        print(command)
        
        #do this via a "bat" file because of diamond bug that does not want to work with custom matrices
        batfile=str(Path(basedir,"alignment.bat"))
        with open(batfile,"w") as bat:
            bat.write("#!/bin/bash"+"\n"+command)

        stdout, stderr =subprocess.Popen(batfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        shutil.move(str(Path(basedir,"diamond.log")), str(Path(Output_directory ,output_folder,Path(input_file).stem+".log")))
        output_files.extend([output_file])
        
    return output_files  
