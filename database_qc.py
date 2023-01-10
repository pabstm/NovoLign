# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:08:28 2022
@author: ZR48SA
modified MP 12 Dec 2022

"""
from __main__ import *
from process_alignment import lca
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

def assign_ncbi_taxonomy(database_psm_file,# df containing columns: Scan, Accessions (with Scan being a list of Accessions matched to that scan)
                        fasta_file,        # path to fasta file used for annotation
                        method="OX",       # "OX", based on taxon id, or "OS", based on name                       
                        acc_dliml="sp|",   # left delimiter of accession information in fasta header
                        acc_dlimr=" ",     # right delimiter of accession information in fasta header                       
                        tax_dliml="OX=",   # left delimiter of taxonomic information in fasta header
                        tax_dlimr=" ",     # right delimiter of taxonomic information in fasta header  
                        ):
    
    db_df=pd.read_csv(database_psm_file)
    db_df["Accession"]=db_df["Accession"].str.split(":")
    odf=db_df[["Scan","Accession"]]
    if type(odf["Accession"].iloc[0])==list:
        odf=odf.explode("Accession")
    adf=[]
    for record in SeqIO.parse(fasta_file,format="fasta"):
        d=record.description
        adf.append([d.split(tax_dliml)[1].split(tax_dlimr)[0],
                   d.split(acc_dliml)[1].split(acc_dlimr)[0]])
    adf=pd.DataFrame(adf,columns=[method,"Accession"])
    odf=odf.merge(adf,on="Accession")
    odf=odf.merge(ncbi_taxdf.astype(str),on=method,)
    ldf=lca(odf[["Scan"]+ranks],group_on="Scan") # if single lca return
    out=db_df.merge(ldf,on="Scan",how="left")
    return out 


def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0
    for n in range(0, N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n
    result *= x ** (N + 1)
    return result


class rectangle:
    def __init__(self,taxa,count,color):
        self.taxa=taxa
        self.count=count
        self.color=color # good/dump/gap/missing


def merge_taxonomy(dn_df,
                   db_df):
    dn=dn_df[["Scan"]+ranks].astype(str)
    dn.columns=["Scan"]+["dn_"+rank for rank in ranks]
    db=db_df[["Scan"]+ranks].astype(str)
    db.columns=["Scan"]+["db_"+rank for rank in ranks]
    merged_taxonomy=dn.merge(db,on="Scan",how="outer").fillna("nan")
    return merged_taxonomy


def Topx_taxa(merged_taxonomy,rank,topx=10):
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
    

def Compare_Bar(denovo_peptides_lca,database_searching_file,fasta_database,
                target_ranks=["genus"],
                write_figure=True,
                write_data=True):
    
    database_peptides=assign_ncbi_taxonomy(database_searching_file,fasta_database)
    merged_taxonomy=merge_taxonomy(denovo_peptides_lca,database_peptides)

    target_ranks=["genus"]
    cols=["dn_"+rank for rank in ranks]+["db_"+rank for rank in ranks]

    for rank in target_ranks:        
        # get topx taxa
        unitax,colordict=Topx_taxa(merged_taxonomy,rank)
        
        # make stacked bar graph
        jacols=["dn_"+rank,"db_"+rank]
        ois=[]
        for jcol in jacols:
            oi=merged_taxonomy[jcol]
            oi=oi[oi.isin(unitax)].value_counts()
            ois.append(oi)
        ois=pd.concat(ois,axis=1).fillna(0)
        
        # get DN only
        dn_only=merged_taxonomy.loc[merged_taxonomy.db_superkingdom=="nan",["dn_"+rank]]
        dn_only=dn_only[dn_only["dn_"+rank].isin(unitax)].value_counts().to_frame().reset_index().set_index("dn_"+rank)
        # get LQ DB only
        dblq_only=merged_taxonomy.loc[merged_taxonomy.dn_superkingdom=="nan",["db_"+rank]]
        dblq_only=dblq_only[dblq_only["db_"+rank].isin(unitax)].value_counts().to_frame().reset_index().set_index("db_"+rank) 
        # make bars
        ois=pd.concat([dn_only,ois,dblq_only],axis=1).fillna(0)
        ois.columns=["DN_only","DN_all","DB_all","DBLQ_only"]
        nois=ois/ois.sum()*100
        # make path
        output_folder=str(Path(Output_directory,"database_qc")) 
        if not os.path.exists(output_folder): os.mkdir(output_folder)
        
        # create graph
        titles=["absolute","normalized"]
        ylabels=["# of scans","norm abundance"]
        for ix,i in enumerate([ois,nois]):
            n=len(i)
            cmap=sns.color_palette("hls", n_colors=n)
            fig=i.T.plot(kind='bar', stacked='True', width=.9,color=cmap)
            plt.ylabel(ylabels[ix])
            plt.title(titles[ix])
            sns.move_legend(fig, "upper left", bbox_to_anchor=(1,1))    

            if write_figure:
               fig1=fig.get_figure()
               fig1.savefig(str(Path(output_folder,rank+"_"+Path(de_novo_file).stem+"_"+str(titles[ix])+"_topX.png")),dpi=400,bbox_inches="tight")
               
        if write_data:
            ois.to_csv(str(Path(output_folder,rank+"_"+Path(de_novo_file).stem+"_topx_bars.tsv")),sep="\t")
            merged_taxonomy.to_csv(str(Path(output_folder,rank+"_"+Path(de_novo_file).stem+"_DN_DB_merged.tsv")),sep="\t")
