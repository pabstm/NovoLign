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
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()

#%% import 
import urllib, gzip, zipfile, shutil, tarfile, time, requests
from sys import platform


#%% Download single file

def download_extract(urls,path):
    
    if type(urls)==str:
        urls=[urls]
    
    for url in urls:
        print(url)
        
        if not os.path.exists(path): os.mkdir(path)
        filename  = str(Path(path,url.split("/")[-1]))
        
        #download
        if not os.path.exists(filename): #check if file is already downloaded
            while True:
                try:
                    urllib.request.urlretrieve(url,filename)
                    break
                except:
                    time.sleep(2)
                    print("retry")
        
        if not os.path.exists(Path(filename).stem): #check if file is already there extracted
            #recursive extraction            
            while any([f.endswith((".zip",".gz",".tar")) for f in os.listdir(path)]):
                
                for f in os.listdir(path):
                    
                    i=str(Path(path,Path(f)))
                    o=str(Path(path,Path(f).stem))
                    
                    if f[0].isalnum() and f.endswith(".zip"):
                        print("extracting "+f)
                        with zipfile.ZipFile(i, 'r') as zip_ref:
                            zip_ref.extractall(path)
                        if os.path.exists(i): os.remove(i)
        
                    if f[0].isalnum() and f.endswith(".gz"):
                        print("extracting "+f)
                        with gzip.open(i,'rb') as f_in:
                            with open(o,'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                        if os.path.exists(i): os.remove(i)
                        
                    if f[0].isalnum() and f.endswith(".tar"):
                        print("extracting "+f)
                        tar = tarfile.open(i, "r:")
                        tar.extractall(path)
                        tar.close()
                        if os.path.exists(i): os.remove(i)
                


#Wrapper for downloading db
def download_db(DB,path=False):

    #RefSeq
    if DB=="RefSeq":
        files=requests.get("https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein-prot-metadata.json").json().get("files")
        if not path: path=str(Path(basedir,"RefSeq"))
        for file in files:
            download_extract(file,path)
    
    #NCBI NR
    if DB=="NCBI_NR":
        files=requests.get("https://ftp.ncbi.nlm.nih.gov/blast/db/nr-prot-metadata.json").json().get("files")
        if not path: path=str(Path(basedir,"NCBI_NR"))
        for file in files:
            download_extract(file,path)
    
    #Swiss-Prot
    if DB=="Swiss-Prot":
        url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        if not path: path=str(Path(basedir,"Swiss-Prot"))
        download_extract(url,path)
    
    #TrEMBL
    if DB=="TrEMBL":
        url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
        if not path: path=str(Path(basedir,"TrEMBL"))
        download_extract(url,path)
    
    #UniProt (Swiss-Prot+TrEBML)
    if DB=="UniProt":
        url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        if not path: path=str(Path(basedir,"UniProt"))
        download_extract(url,path)
        
        url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
        if not path: path=str(Path(basedir,"UniProt"))
        download_extract(url,path)
    
    #Uniref 100
    if DB=="UniRef100":
        url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
        if not path: path=str(Path(basedir,"UniRef100"))
        download_extract(url,path)
    
    #Uniref 90
    if DB=="UniRef90":
        url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
        if not path: path=str(Path(basedir,"UniRef90"))
        download_extract(url,path)
    
    #Uniref 50
    if DB=="UniRef50":
        url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
        if not path: path=str(Path(basedir,"UniRef50"))
        download_extract(url,path)
    
    
    return path



def download_ncbi_taxdump(path=False):

    url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"    
    if not path: path=str(Path(basedir,"taxdump"))
    download_extract(url,path)
    
    #cleanup
    for i in os.listdir(path):
        if i=="names.dmp":
            names=str(Path(path,i))
        elif i=="nodes.dmp":
            nodes=str(Path(path,i))
        else:
            os.remove(i)

    return path,names,nodes


def download_diamond(path=False):


    if platform == "linux" or platform == "linux2":
        # linux
        url="https://github.com/bbuchfink/diamond/releases/download/v2.0.12/diamond-linux64.tar.gz"
    elif platform == "win32":
        # Windows
        url="https://github.com/bbuchfink/diamond/releases/download/v2.0.12/diamond-windows.zip"    
    else:
        print("this pipeline only supports widows or linux systems, since DIAMOND installation on OSX requires Miniconda")
        return
    
    if not path: path=str(Path(basedir,"diamond"))
    download_extract(url)
    
    for i in os.listdir(path):
        if i.starswith("diamond"):
            return str(Path(path,i))
    
    

