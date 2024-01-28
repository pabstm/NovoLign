# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:08:28 2022
@author: ZR48SA
modified MP 12 Dec 2022

"""
ranks=["superkingdom","phylum","class","order","family","genus","species"]
from process_alignment import lca
import matplotlib.pyplot as plt
from collections import Counter
import more_itertools as mit
from pathlib import Path
from Bio import SeqIO
import seaborn as sns
import pandas as pd
import numpy as np
import os

def assign_ncbi_taxonomy(database_psm_file,# df containing columns: Scan, Accessions (with Scan being a list of Accessions matched to that scan)
                        fasta_file,        # path to fasta file used for annotation
                        ncbi_taxdf,        
                        taxids='',         # Diamond taxid output file
                        method="OX",       # "OX", based on taxon id, or "OS", based on name                       
                        acc_dlimr=" ",     # Default right delimiter of accession information in fasta header                       
                        #acc_dlimr="|",    # A baumanii, Aeromonas, Strep mutans, Kleiner, WWTP UP, Chl, HCR, Para ...
                        tax_dliml="OX=",   # left delimiter of taxonomic information in fasta header
                        tax_dlimr=" ",     # right delimiter of taxonomic information in fasta header  
                        ):
    
        db_df=pd.read_csv(database_psm_file)
        db_df["Accession"]=db_df["Accession"].str.split(":")
        odf=db_df[["Scan","Accession"]]

        if type(odf["Accession"].iloc[0])==list:
            odf=odf.explode("Accession")
        odf.Accession=odf.Accession.str.replace('sp\|', "").str.replace('tr\|', "")
            
        adf=[]
        if len(taxids)<1: #use fasta "OX" info
            for record in SeqIO.parse(fasta_file,format="fasta"):
                d=record.description
                adf.append([d.split(tax_dliml)[1].split(tax_dlimr)[0], d.split('sp|')[1].split(acc_dlimr).pop(0) if 'sp|' in d else d.split('tr|')[1].split(acc_dlimr).pop(0)])
            adf=pd.DataFrame(adf,columns=[method,"Accession"])
            odf=odf.merge(adf,on="Accession")
            
        else: #use Diamond taxid file from in input folder
            taxids.Accession=taxids.Accession.str.replace('sp\|', "").str.replace('tr\|', "")
            #taxids["Accession"] = taxids["Accession"].str.split('|').str[0] # used for P den
            odf=odf.merge(taxids,on="Accession")
            
        odf=odf.merge(ncbi_taxdf.astype(str),on=method,)
        ldf=lca(odf[["Scan"]+ranks],Output_directory="simple_lca",denovo_peptides="simple_lca",group_on="Scan",method="standard")
        out=db_df.merge(ldf,on="Scan",how="left")
        return out


class rectangle:
    def __init__(self,taxa,count,color):
        self.taxa=taxa
        self.count=count
        self.color=color # good/dump/gap/missing


def merge_taxonomy(dn_df,
                   db_df):
    dn=dn_df[["Scan"]+ranks].astype(str)
    dn.columns=["Scan"]+["dn_"+rank for rank in ranks]
    db_df['psm']='Y' #add db search identifier
    db=db_df[["psm"]+["Scan"]+ranks].astype(str)
    db.columns=["psm"]+["Scan"]+["db_"+rank for rank in ranks]
    merged_taxonomy=db.merge(dn,on="Scan",how="outer").fillna("nan")
    return merged_taxonomy


def Topx_taxa(merged_taxonomy,rank,topx=15):
        xdf=pd.concat([
        pd.DataFrame(Counter(merged_taxonomy["db_"+rank]).most_common(),columns=[rank,"db"]).set_index(rank),
        pd.DataFrame(Counter(merged_taxonomy["dn_"+rank]).most_common(),columns=[rank,"dn"]).set_index(rank)],
            axis=1)
        xdf=xdf/xdf.sum()
        uxdf=xdf.loc[[i for i in xdf.index if (i!="nan") & (i!="")]]
        unitax=uxdf.sum(axis=1).sort_values()[::-1][0:topx].index # select topx most frequent taxa
        colordict={} 
        cmap=sns.color_palette("Paired",n_colors=topx) 
        [colordict.update({i:cmap[c]}) for c,i in enumerate(unitax)]  
        return unitax,colordict


def fill_g(x):
    cons=[np.array(list(g)) for g in mit.consecutive_groups(np.argwhere(x==""))]
    if cons:
        for c in cons:  
            c=c.flatten().tolist()
            if c[-1]!=len(x)-1:
                for i in c:
                    x[i]="gap_"+str(x[c[-1]+1])+"_"+ranks[i]
    return list(x)


def Compare_Bar(denovo_peptides_lca,database_searching_file,fasta_database,taxa,Output_directory,ncbi_taxdf,
                target_ranks=["order","family","genus"],
                write_figure=True,
                write_data=True,
                fillgaps=True
                ):
    #%%
    
    
    dn_lca=pd.read_csv(denovo_peptides_lca, delimiter='\t')
    
    if taxa=='':
        database_peptides=assign_ncbi_taxonomy(database_searching_file,fasta_database,ncbi_taxdf)
    else:
        taxids = pd.read_csv(taxa, sep = '\t', header=None)
        taxids.columns = ["Accession","OX", "e-value"]
        taxids['OX'] = taxids['OX'].apply(str)
        database_peptides=assign_ncbi_taxonomy(database_searching_file,fasta_database,ncbi_taxdf,taxids=taxids)
        database_peptides=database_peptides.fillna("")
    
    if fillgaps==True:
        gaps=database_peptides[(database_peptides[ranks]=="").any(axis=1)]
        u,ix,inv=np.unique(gaps[ranks].apply(";".join,axis=1),return_inverse=True,return_index=True)
        u=gaps.iloc[ix][ranks].values
        gap_1=pd.DataFrame(np.array(list((map(fill_g,u))))[inv]).set_index(gaps.index)
        if gaps.empty == False:
            database_peptides.loc[gap_1.index,ranks]=gap_1.values      
            
    merged_taxonomy=merge_taxonomy(dn_lca,database_peptides)


    #make bar graph
    for rank in target_ranks:        
        unitax,colordict=Topx_taxa(merged_taxonomy,rank)
        jacols=["dn_"+rank,"db_"+rank]
        ois=[]
        for jcol in jacols:
            oi=merged_taxonomy[jcol]
            oi=oi[oi.isin(unitax)].value_counts()
            ois.append(oi)
        ois=pd.concat(ois,axis=1).fillna(0)
        
        # get DN only
        dn_only=merged_taxonomy.loc[merged_taxonomy.psm=="nan",["dn_"+rank]]
        dn_only=dn_only[dn_only["dn_"+rank].isin(unitax)].value_counts().to_frame().reset_index().set_index("dn_"+rank)
        # get LQ DB only
        dblq_only=merged_taxonomy.loc[merged_taxonomy.dn_superkingdom=="nan",["db_"+rank]]
        dblq_only=dblq_only[dblq_only["db_"+rank].isin(unitax)].value_counts().to_frame().reset_index().set_index("db_"+rank) 
        # make bars
        ois=pd.concat([dn_only,ois,dblq_only],axis=1).fillna(0)
        ois.columns=["DN_only","DN_all","DB_all","DB_only"]
        ois=ois[["DB_all","DB_only","DN_all","DN_only"]]
        
        nois=ois/ois.sum()*100
        # make path
        output_folder=str(Path(Output_directory,"database_qc")) 
        if not os.path.exists(output_folder): os.mkdir(output_folder)
        
        # create graph
        titles=["absolute","normalized"]
        ylabels=["Number of scans","Normalized abundance"]
        for ix,i in enumerate([ois,nois]):
            n=len(i)
            cmap=sns.color_palette("hls", n_colors=n)
            fig=i.T.plot(kind='bar', stacked='True', width=.9,color=cmap)
            plt.ylabel(ylabels[ix])
            plt.title(titles[ix])
            sns.move_legend(fig, "upper left", bbox_to_anchor=(1,1))    

            if write_figure:
               fig1=fig.get_figure()
               fig1.savefig(str(Path(output_folder,"DB_vs_DN_"+rank+str(titles[ix])+"_topX.png")),dpi=400,bbox_inches="tight")
               
            if write_data:
                ois.to_csv(str(Path(output_folder,"DB_vs_DN_"+rank+"_topx_bars.tsv")),sep="\t")
    if write_data:
        merged_taxonomy.to_csv(str(Path(output_folder,"DN_DB_merged.tsv")),sep="\t")
