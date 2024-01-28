"""
module for taxonomic grouping using various LCAs

Author: mpabst
Date: 26-06-2022
"""

from pathlib import Path
import pandas as pd
import numpy as np
from itertools import chain
from collections import Counter
import glob, os

def stacked_bar(ranks,df,ylabel,pathout,filename,plt,os,np): 

    labels=[i.split("_name")[0] for i in ranks];
    countcols=[i for i in df.columns if "count" in i]
    absdat=df[countcols] 
    absdat.columns=labels
    normdat= absdat/np.nansum(absdat.to_numpy(),axis=0)
    
    figure, axes = plt.subplots(1, 2)
    ax1=absdat.T.plot(ax=axes[0],kind='bar', stacked=True, figsize=(10, 6), legend=False)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel('taxonomic ranks')
    ax1.set_xticklabels(labels, rotation=30)
    ax1.set_title("Absolute")
    
    ax2=normdat.T.plot(ax=axes[1],kind='bar', stacked=True, figsize=(10, 6), legend=False)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel('taxonomic ranks')
    ax2.set_xticklabels(labels, rotation=30)
    ax2.set_title("Normalized")

    plt.gcf().suptitle(Path("test").stem) # placeholder basename
    figname=str(Path(pathout,(filename.replace(os.path.splitext(filename)[1], '.png'))))
    plt.savefig(figname)
    
    return figname

def graphs(file,Output_directory,ranks,cutbranch):
    
    tsv=pd.read_csv(file, sep='\t', header=0)
    denoise="off" # switched off because postfilter used after lca
    cutbranch=5 # default cutoff frequency
    tax_normalize=False # normalize to total for that rank
    xlsdf=[]
    xlsdf=tsv
    
    # tax frequency filter    
    if denoise=="on":
        unirows=xlsdf[["Peptide"]+ranks].drop_duplicates()[ranks].astype(str).values.tolist()
        jbranch=["#*".join(i) for i in unirows]
        fbranch=[branch for branch, counts in Counter(jbranch).items() if counts >= cutbranch]
        allowed_taxids=set(chain(*[i.split("#*") for i in fbranch]))
        for i in ranks:
            xlsdf.loc[~xlsdf[i].isin(allowed_taxids),i]="" 

    # create output folder
    pathout=str(Path(Output_directory,"composition"))
    if not os.path.exists(pathout): os.makedirs(pathout)     

    # tax composition
    quantdf=pd.DataFrame()
    for rank in ranks:
       
        values=Counter(xlsdf[rank].astype(str))
        if "" in values.keys(): values.pop("")
            
        values=pd.Series(values).sort_values(axis=0, ascending=False, inplace=False).reset_index()
        values.columns=[rank,rank+"_count"]   
        values=values[values[rank]!=""]
        if tax_normalize==True: 
            values[rank+"_count"]=values[rank+"_count"]/values[rank+"_count"].sum()*100
        
        quantdf=pd.concat([quantdf, values], axis=1)
        
    # get name extension
    ext=Output_directory.split("Output_")[1]
    
    # write composition
    quantdfs=quantdf.fillna(0)
    namecols=[i for i in quantdf.columns  if "count" not in i] 
    countcols=[i for i in quantdf.columns  if "count" in i] 
    name_1=file.split('tab\\')
    basename="comp_Cts_"+Path(name_1[0]).stem+'_'+ext+'.xlsx'
    xlsfilename=str(Path(pathout,basename))
    writer=pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
    quantdf[namecols].to_excel(writer, sheet_name='TAX_LINEAGES')
    quantdf[countcols].to_excel(writer, sheet_name='TAX_COUNTS')
    quantdf.to_excel(writer, sheet_name='COMBINED')  
    
    #depends on pandas version
    try:
        writer.save() 
    except Exception as error:
        writer.close()

    return quantdfs
