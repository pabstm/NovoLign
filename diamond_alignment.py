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
                      matrix_path=str(Path(basedir,"PAM70.mat")),
                      output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"],
                      select=" -k50 ", # can be --top x% or --max-target-seqs/-k
                      block_size=5,
                      index_chunks=1,
                      minimum_pident=80,
                      minimum_coverage=80,
                      minimum_bitscore=20,
                      gap_open=0,
                      gap_extend=8 
                      ):
    


    output_folder=str(Path(Output_directory,"diamond_alignments"))    
    if not os.path.exists(output_folder): os.mkdir(output_folder)
    
    output_files=[]
    for input_file in input_files:
        output_file=str(Path(output_folder,Path(input_file).stem+".tsv"))
        
        command="cd "+'"'+basedir +'"'+ " && " + \
                "".join(['"'+str(Path(basedir,"diamond"))+'"',   
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
        
        # print
        stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        print(stderr) 
        shutil.move(str(Path(basedir,"diamond.log")), str(Path(Output_directory,output_folder,Path(input_file).stem+".log")))
        
        output_files.extend([output_file])
        
    return output_files  
