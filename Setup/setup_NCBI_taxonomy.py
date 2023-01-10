
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:34:40 2021

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

import os
from pathlib import Path
import pandas as pd
import numpy as np





#ranks to be included

#%% iteratively construct taxonomies from nodes.dmp

def parse_NCBI_taxonomy(names,nodes,ranks=["superkingdom","phylum","class","order","family","genus","species"],path=False):

    if not path: path=str(Path(basedir,"ncbi_taxonomy"))    
    if not os.path.exists(path): os.mkdir(path)

    with open(nodes,"r") as f: lines=pd.DataFrame([l.split("\t")[0:5] for l in f.readlines()])
    nodedf=lines.iloc[:,[0,2,4]]
    nodedf.columns=["taxid","parent_taxid","rank"]
    
    
    #has_rank=nodedf[nodedf["rank"]==ranks[-1]] #select those that have the final rank (in this case species)
    has_rank=nodedf
    
    
    
    inc=0               #used for renaming columns
    count=len(has_rank) #used for breaking option 1 (no more remaining candidates)
    break_counter=10     #used for breaking option 2 (3 loops same output)
    
    
    parent=nodedf.copy()
    
    completed=[]
    while count!=0:
        
        
    
        #renaming
        parent.columns=[str(inc)+"_taxid",str(inc)+"_parent_taxid",str(inc)+"_rank",]
        
        if inc:
            has_rank=has_rank.merge(parent,left_on=str(inc-1)+"_parent_taxid",
                                    right_on=str(inc)+"_taxid",how="left")
            
        else:
            has_rank=has_rank.merge(parent,left_on="parent_taxid",
                                    right_on=str(inc)+"_taxid",how="left")
            
    
    
    
        d=(has_rank[str(inc)+"_rank"]==ranks[0]).sum() #to check progress, substract those that have the first rank
        
        ix=has_rank[has_rank[str(inc)+"_rank"]=="superkingdom"].index
        
        
        completed.append(has_rank.loc[ix.tolist(),:])
        has_rank=has_rank.drop(ix)
        #this is to break out of the loop if there are no new hits within x iterations
        if d==0:
            break_counter+=1
        else:
            break_counter=0
        if break_counter>=10:
            break
        
    
    
        count-=d
        inc+=1
        print(count)
        
    
    rsp=[]
    for ix,has_rank in enumerate(completed):
        print("parsing batch "+str(ix)+"/"+str(len(completed)-1))
        rs=[]
        for r in ranks:
            print(r)
            v=np.argwhere((has_rank==r).values)
            has_rank.values[v[:,0],v[:,1]]
            rdf=pd.DataFrame([has_rank.values[v[:,0],0], has_rank.values[v[:,0],v[:,1]-2]]).T
            rdf.columns=["idx",r]
            rdf=rdf.set_index("idx")
            rs.append(rdf)
        
        p=rs[0]
        for r in rs[1:]:
            p=p.merge(r,how="left",on="idx")
        rsp.append(p)
    
    p=pd.concat(rsp)
    
        
    #rename taxids according to names.damp
    
    with open(names,"r") as f: lines=pd.DataFrame([l.split("\t")[0:7] for l in f.readlines()])
    namesdf=lines.iloc[:,[0,2,4,6]]
    namesdf.columns=["taxid","name","x","type"]
    namesdf=namesdf.loc[namesdf["type"]=="scientific name",["taxid","name"]]
    namesdf=namesdf.set_index("taxid")
    
    #root is used as placeholder for nans
    oi=pd.DataFrame(namesdf.loc[p.fillna("1").values.flatten()].values.reshape(-1,len(ranks)),columns=ranks) 
    oi[oi=="root"]=""
    oi.index=p.index
    #write outputs
    p=p.fillna("")
    #p.to_csv("taxids.tsv",sep="\t")
    
    oi.index.name="OX"
    namesdf.index.name="OX"
    namesdf.columns=["OS"]
    oi=oi.merge(namesdf,on="OX")
    
    outfile=str(Path(path,"parsed_ncbi_taxonomy.tsv"))
    oi.to_csv(outfile ,sep="\t")
    
    return outfile

