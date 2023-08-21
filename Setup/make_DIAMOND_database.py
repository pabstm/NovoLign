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


import subprocess
def make_diamond_database(diamond_path,database_path,outpath,path=False,delete_old=False):
    
    # if not path:
    #     path=database_path
    # ppath=Path(path).parents[0]
    # #outpath=str(Path(ppath,Path(ppath.stem))) #remove path extension (needed for diamond command syntax)

    command="cd "+'"'+str(Path(diamond_path).parents[0]) +'"'+ " && "
    command+="diamond makedb --in "+'"'+database_path+'"' + " -d "+'"'+outpath+'"'
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    if delete_old:
        os.remove(database_path)

    return outpath+".dmnd"


