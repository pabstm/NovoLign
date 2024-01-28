# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:33:27 2022
@author: ZR48SA
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
import seaborn as sns
import pandas as pd
import numpy as np

def Plot_high_scoring(peptide_df, # all possible DN sequences
                      target, # all aligned target DN sequences
                      database_searching_file, # DB searching psms
                      Output_directory,
                      de_novo_file,
                      metric=None,
                      write_figure=True,
                      write_data=True,
                      ):
        
    output_folder=str(Path(Output_directory,"experiment_qc")) 
    if not os.path.exists(output_folder): os.mkdir(output_folder)
    target=target.sort_values(by=["Peptide","score"])
    best=target[["Peptide","qseq","sseq","score","pident"]].groupby("Peptide").nth(0).reset_index()
     
    if metric==None:
       metric=[i for i in ["predicted_score",'ALC (%)'] if i in peptide_df.columns][0] #which collumn to plot
      
    best["Seq_len"]=best["Peptide"].apply(len)
    best["q_len"]=best.qseq.apply(len)
    best["s_len"]=best.sseq.apply(len)
    best["query_cover"]=best["s_len"]/best["Seq_len"]
     
    best["Match"]="aligned"
    best.loc[best["query_cover"]<1,"Match"]="aligned tag"
    best.loc[best["qseq"]==best["sseq"],"Match"]="exact tag"
    best.loc[best["Peptide"]==best["sseq"],"Match"]="exact"
    m=peptide_df.merge(best[["Peptide","Match"]],on="Peptide",how="left")
    m[metric]=m[metric].astype(float)
    m=m.fillna("unmatched")
     
    # make bins
    bw=5
    bins=np.arange(m[metric].min(),m[metric].max(),bw)
    m["bin"]=np.digitize(m[metric],bins)
    g=m.groupby(["bin","Match"]).size()
    pv=pd.pivot_table(g.reset_index(),values=[0], index=['bin'],columns=['Match'])
    pv.columns=pv.columns.droplevel(0)
    pvn=pv.divide(pv.sum(axis=1),axis=0)
    pvn.index=bins
    cols=[i for i in ["exact","exact tag","aligned","aligned tag","unmatched"] if i in pvn.columns]
    pvn=pvn[cols]
    pvn[cols[:-1]].sum()
    pvn["fraction_matched"]=pvn[cols[:-1]].sum(axis=1)/(pvn[cols[:-1]].sum(axis=1)+pvn["unmatched"])
     
    # plot alignment histogram
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(1,1,1)
     
    hplot=sns.histplot(data=m, ax=ax, stat="count", multiple="stack",
                 common_norm=True,
                 x=metric, kde=False,
                 hue="Match",
                 binwidth=bw,
                 palette=["#9FA5A9","#992714","#D5671E","#D5C31E","#8CD51E"],
                 hue_order=cols[::-1],
                 element="bars", legend=True)
    ax.set_title(Path(de_novo_file).stem)
    ax.set_xlabel(metric)
    ax.set_ylabel("Number of scans")

    # plot database searching histogram
    if database_searching_file!=None:
        db_df=pd.read_csv(database_searching_file)
    
        db_m=m.copy()
        db_m["Match"]="not detected"
        db_m.loc[db_m.Scan.isin(db_df.Scan),"Match"]="Matched"
        db_m=db_m[["Match",metric]]
        
        figdb = plt.figure(figsize=(7,5))
        ax = figdb.add_subplot(1,1,1)
         
        hplot=sns.histplot(data=db_m, ax=ax, stat="count", multiple="stack",
                      common_norm=True,
                      x=metric, kde=False,
                      hue="Match",
                      binwidth=bw,
                      palette=["#9FA5A9","#992714"],
                      hue_order=["Matched","not detected"][::-1],
                      element="bars", legend=True)
        ax.set_title(Path(database_searching_file).stem)
        ax.set_xlabel(metric)
        ax.set_ylabel("Number of scans")
    
        if write_figure:
            figdb.savefig(str(Path(Output_directory,output_folder,Path(database_searching_file).stem+"_bins.png")),dpi=300,bbox_inches="tight")
                                                                  
    # mass error alc scatter plot
    cmap = plt.get_cmap("tab10")
    category=["exact","exact tag","aligned","aligned tag","unmatched"][::-1]
    palette=["#9FA5A9","#992714","#D5671E","#D5C31E","#8CD51E"]

    fig1 = plt.figure(figsize=(7, 6))
    fig1.suptitle(Path(de_novo_file).stem)
    grid = plt.GridSpec(4, 4, hspace=0.3, wspace=0.4)
    main_ax = fig1.add_subplot(grid[:-1, 1:])
    y_hist = fig1.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
    x_hist = fig1.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)

    for ix,cat in enumerate(category):
        d=m[m["Match"]==cat]
        main_ax.plot(d["ALC (%)"].astype(float),d["ppm"],'.',markersize=0.4,color=palette[ix],alpha=0.8)
        yall,x,_=y_hist.hist(d["ppm"],bins=50, density=True, histtype='stepfilled', orientation='horizontal',color=palette[ix], alpha=0.4)
        yall,x,_ =x_hist.hist(d[metric].astype(float), bins=50, density=True,histtype='stepfilled',  orientation='vertical', color=palette[ix], alpha=0.4)

    y_hist.set_ylabel("ppm")
    y_hist.set_yticks(np.arange(m["ppm"].min(),m["ppm"].max(),5))
    y_hist.set_yticklabels(y_hist.get_yticks(),fontsize=8)

    x_hist.set_xlabel(metric)
    x_hist.set_xticks(np.arange(m[metric].min(),m[metric].max(),5))
    x_hist.set_xticklabels(x_hist.get_xticks(),fontsize=8)

    plt.legend(category,loc=[-0.4,-0.1])
    
    if write_figure:
        fig.savefig(str(Path(Output_directory,output_folder,Path(de_novo_file).stem+"_bins.png")),dpi=300,bbox_inches="tight")
        fig1.savefig(str(Path(Output_directory,output_folder,Path(de_novo_file).stem+"_scatter.png")),dpi=300,bbox_inches="tight")
    if write_data:
        pv.to_csv(str(Path(Output_directory,output_folder,Path(de_novo_file).stem+"_matched_bins_absolute.tsv")),sep="\t")
        pvn.to_csv(str(Path(Output_directory,output_folder,Path(de_novo_file).stem+"_matched_bins_normalized.tsv")),sep="\t")
        m.to_csv(str(Path(Output_directory,output_folder,Path(de_novo_file).stem+"_dn_matched_peptides.tsv")),sep="\t")
        db_m.to_csv(str(Path(Output_directory,output_folder,Path(de_novo_file).stem+"_db_matched_peptides.tsv")),sep="\t")