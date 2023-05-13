#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 09:45:27 2022
@author: hugokleikamp
"""
#%% import modules, set paths
from pathlib import Path
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import more_itertools as mit
from itertools import chain
from collections import Counter
from __main__ import ranks, ncbi_taxdf#, de_novo_file
import re

# read alignment file to generator and filter on topx % scoring


def filter_chunk(s,top_score_fraction,minimum_bitscore):
    
    s=s.assign(Peptide=s.qseqid.apply(lambda x: x.split(";")[0]))
    s.loc[:,"qlen"]= s["Peptide"].apply(len)
    s.loc[:,"slen"]= s["sseq"].apply(len)
    s.loc[:,"score"]=s["bitscore"]*s["slen"]/s["qlen"]*s["pident"]/100 #bitscore corrected for alignment length
    
    if top_score_fraction:    s=s[s.score>=s.score.max()*(top_score_fraction)]
    if minimum_bitscore:      s=s[s.bitscore>=minimum_bitscore]
    
    return s


def Diamond_alignment_Reader(input_file,
                             
                             Output_directory,
                             top_score_fraction=0.9,
                             minimum_bitscore=0, 
                             Read_batch=1000000,
                             output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]):
    
    cdf=pd.read_csv(input_file, sep='\t', chunksize=Read_batch,names=output_columns) # read to generator
 
    sc=[]
    print("reading alignment chunks:")
    for ix,c in enumerate(cdf):
        print(ix)
        c["bitscore"]=c["bitscore"].astype(float)
        _,index=np.unique(c.qseqid,return_index=True)
        index.sort()
        index_m=index.tolist()+[max(index)+1]
       
        if ix>0:
            c=pd.concat([sc[-1],c])
       
        sc=[c.iloc[index_m[n]:index_m[n+1]] for n in range(len(index_m)-1)]
        for s in sc[:-1]:
            yield filter_chunk(s,top_score_fraction=top_score_fraction,minimum_bitscore=minimum_bitscore)
        
    # last one
    s=sc[-1]
    yield filter_chunk(s,top_score_fraction=top_score_fraction,minimum_bitscore=minimum_bitscore)
 

def fill_g(x):
    cons=[np.array(list(g)) for g in mit.consecutive_groups(np.argwhere(x==""))]
    if cons:
        for c in cons:  
            c=c.flatten().tolist()
            if c[-1]!=len(x)-1:
                for i in c:
                    x[i]="gap_"+str(x[c[-1]+1])+"_"+ranks[i]
    return list(x)


def Process_diamond_alignment(alignments,
                              
 
                              
                              Output_directory,
                              fill_gaps=True,

                              
                              top_score_fraction=0.9,
                              minimum_bitscore=0, 
                              
                              Taxid_delimiter="OX=",Gene_delimiter="GN=",
                              weight_rank="genus",
                              
                              output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]
                              ):
    
    # #%%
    # #test
    # alignments
    # Output_directory
    # fill_gaps=True

    
    # top_score_fraction=0.9
    # minimum_bitscore=0 
    
    # Taxid_delimiter="OX="
    # Gene_delimiter="GN="
    # weight_rank="genus"
    
    
    target_decoys=[]
    if type(alignments)==str: alignments=[alignments]
    for alignment in alignments:
        

        iterable=Diamond_alignment_Reader(alignment,Output_directory,top_score_fraction=top_score_fraction,minimum_bitscore=minimum_bitscore,output_columns=output_columns)
        al=pd.concat([i for i in iterable])
       
        # add NCBI taxonomy UniRef or UniprotKB headers
        if  al.head(10)["sseqid"].str.contains("Taxid=").any(): Taxid_delimiter,Gene_delimiter="Taxid=","RepID="
        al["OX"]  =al["stitle"].apply(lambda x: ((x+Taxid_delimiter+" ").split(Taxid_delimiter)[1].split()+[" "])[0]).astype(int)
        al["gene"]=al["stitle"].apply(lambda x: ((x+ Gene_delimiter+" ").split( Gene_delimiter)[1].split()+[" "])[0])
        
        al=al.merge(ncbi_taxdf,how="left",on="OX").fillna("")
       
        # filling annotation gaps in lineages
        if fill_gaps:
            gaps=al[(al[ranks]=="").any(axis=1)]
            u,ix,inv=np.unique(gaps[ranks].apply(";".join,axis=1),return_inverse=True,return_index=True)
            u=gaps.iloc[ix][ranks].values
            gap_1=pd.DataFrame(np.array(list((map(fill_g,u))))[inv]).set_index(gaps.index)
            al.loc[gap_1.index,ranks]=gap_1

        al=pd.concat([al, # add information stored in the header of the fasta file to the read dataframe
                  pd.DataFrame.from_dict(al.qseqid.str.split(";").apply(lambda x: eval(x[1])).values.tolist(),orient="columns")],axis=1)

        if "alignment_Target_Decoy" in al.columns:
            target=al[al["alignment_Target_Decoy"]=="Target"]
            decoy=al[al["alignment_Target_Decoy"]=="Decoy"]

            # add weights
            if type(weight_rank)!=type(None):

                # add weights
                mw=pd.concat([d.groupby(weight_rank).size() for d in [target,decoy]],axis=1).fillna(0).reset_index()
                mw.columns=[weight_rank,"target_count","decoy_count"]
                mw["ntarget_count"]=mw["target_count"]/mw["target_count"].sum()
                mw["ndecoy_count"]=mw["decoy_count"]/mw["decoy_count"].sum()
                mw["corrected_weights"]=((mw.ntarget_count-mw.ndecoy_count)/(mw.ntarget_count+mw.ndecoy_count)+1)/2                        # percentile difference normalized to 0-1
                mw["ntarget_count"]=(mw["ntarget_count"]-mw["ntarget_count"].min())/(mw["ntarget_count"].max()-mw["ntarget_count"].min())  # renormalize to 0-1
                mw["ndecoy_count"]=(mw["ndecoy_count"]-mw["ndecoy_count"].min())/(mw["ndecoy_count"].max()-mw["ndecoy_count"].min())       # renormalize to 0-1
                target=target.merge(mw[[weight_rank,"ntarget_count","corrected_weights"]].rename(columns={"ntarget_count":"weights"}),on=weight_rank)
                decoy=  decoy.merge(mw[[weight_rank, "ndecoy_count","corrected_weights"]].rename(columns={ "ndecoy_count":"weights"}),on=weight_rank)

                target["weighted_score"]=target["score"]*target["weights"]
                target["corrected_weighted_score"]=target["score"]*target["corrected_weights"]
                decoy["weighted_score"]=decoy["score"]* decoy["weights"]
                decoy["corrected_weighted_score"]=decoy["score"]* decoy["corrected_weights"]
                
                target_decoys.append(pd.concat([target,decoy]))
         
        else:
            target_decoys.append(al)
        
        
 #%%
    return target_decoys


# Postfilter


def Postfilter(
               lcas,
               denovo_peptides=None,
               group_on="Peptide",      
               filter_cutoff=5,                     # Depending on the input value this has two modes:
                                                    # if >=1, assumes static filter, taxa should appear above this treshold frequency, 
                                                    # if between 0-1 it is used as frequency quantile from the decoy LCA
                ):
 
    if type(denovo_peptides)!=type(None):
        pdf=denovo_peptides[[group_on,'alignment_Target_Decoy']].merge(lcas,how="left",on=group_on)
    else:
        pdf=lcas
    
    unirows=pdf[[group_on]+ranks].drop_duplicates()[ranks].astype(str).values.tolist()
    jbranch=["#*".join(i) for i in unirows]
    
    if filter_cutoff<1 and type(denovo_peptides)!=type(None): #if between 0-1, use decoy quantile
        d=pdf[pdf["alignment_Target_Decoy"]=="Decoy"]
        if len(d): filter_cutoff=d.groupby(ranks).size().quantile(filter_cutoff)
        else: print("warning, no Decoy detected!, select  filter cutoff of 1 or higher")
            
    print("post lca filter: "+str(round(filter_cutoff,2)))
    fbranch=[branch for branch, counts in Counter(jbranch).items() if counts <= filter_cutoff]
    removed_taxons=set(chain(*[i.split("#*") for i in fbranch]))   
    
    for i in ranks:
        lcas.loc[lcas[i].isin(removed_taxons),i]=""
        
    return denovo_peptides,lcas,removed_taxons
    

    
def Taxids_to_fasta(lineages,                 #List of accepted taxonomies taxon names, or tax ids 
                    Database=database_path,   #Database from which to parse the 
                    output_path=db_out,       #output path
                    minimum_rank="OX",        #Minimum rank that supplied taxonomic lineages should have
                    add_decoy=True,           #add reverse decoy of aligned proteins to database
                    add_denovo_peptides=True, #add de novo peptides to database (not included in decoy)
                    
                        ):

        recs=SeqIO.parse(Database,format="fasta")
        chunks=chunk_gen(recs)
        
        if minimum_rank=="OX": #select based directly on supplied taxonomy IDs (simplest approach)
            allowed_taxids=[int(i) for i in lineages if i.isdigit()] 
        
        with open(output_path,"w") as f: #clears file if exists
            pass
        with open(output_path,"a") as f:
            
            for ic,c in enumerate(chunks):
                print("chunk "+ str(ic))
        
                chunk_df=pd.DataFrame([[str(r.seq),r.description,r.id] for r in c],columns=["seq","description","id"])
                taxids=chunk_df["description"].str.split("OX=").apply(lambda x: x[-1].split(" ")[0]).astype(int)

                if minimum_rank!="OX" #merge with ncbi taxonomy and select on a specific rank
                    f=ncbi_taxdf[ncbi_taxdf["OX"].isin(taxids.tolist())]
                    allowed_taxids=f.loc[f[minimum_rank].isin(lineages),"OX"]
                
                chunk_df=chunk_df[taxids.isin(allowed_taxids)]
    
                f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")
                if add_decoy:
                    f.write("\n".join(">decoy_"+chunk_df["description"]+"\n"+chunk_df["seq"].str[::-1)+"\n")
    
        return output_path
    
def Proteins_to_fasta(df,                       #processed diamond alignment 
                      lcas=None,                #input: lca df, used to generate list of proteins to include in database
                      removed_proteins=[],      #input: list of proteins that are excluded from database
                      minimum_rank="",          #input: taxonomic rank, in case lca is supplied, this can only select lcas with a certain level of taxonomic detail 
                      add_decoy=True,           #add reverse decoy of aligned proteins to database
                      add_denovo_peptides=True, #add de novo peptides to database (not included in decoy)
                      ):
    #%%
    if "full_sseq" not in df.columns:
        print("No full protein sequences found, run diamond alignment with full_sseq in output columns!")
        return
   # #%%
   
   #  #test
   #  target_decoys[0]
   #  lcas=None                #input: lca df, used to generate list of proteins to include in database
   #  minimum_rank=""          #input: taxonomic rank, in case lca is supplied, this can only select lcas with a certain level of taxonomic detail 
   #  add_decoy=True           #add reverse decoy of aligned proteins to database
   #  add_denovo_peptides=True
    
    #narrow down database based on lca     
    if type(lcas)!=type(None):
        if minimum_rank in lcas.columns:
            lcas=lcas[lca[minimum_rank]!=""]
        if "proteins" in lcas.columns:
            proteins=list(set(sum(lcas["proteins"].dropna().str.split(", ").tolist(),[])))
            if len(removed_proteins):
                proteins=list(set(proteins) - set(removed_proteins))
            df=df[df.sseqid.isin(proteins)]
    
    
    fasta_str="\n".join((">"+df.stitle+"\n"+df.full_sseq))+"\n"
    if add_decoy:
        fasta_str+="\n".join((">"+"decoy_"+df.stitle+"\n"+df.full_sseq.str[::-1]))+"\n"  #reversed decoy
        
    if add_denovo_peptides:
        dn_peptides=df.qseq.drop_duplicates().reset_index()
        fasta_str+="\n".join((">"+"DeNovo_"+dn_peptides["index"].astype(str)+"\n"+dn_peptides["qseq"]))+"\n"
  
    return fasta_str #long string in fasta format

    

def lca(df,                                  #dataframe with at least a column called Peptide, and ranks
        
        Output_directory,
        group_on="Peptide",                  
        denovo_peptides=None,                #original PSM DataFrame, used to merge back results
        
        #Weighing parameters
        method="weighted",                   # Options: False (conventional LCA), "weighted" (weighted lca), "focused" (focused lca)
        weight_rank="genus",                 # weighted weighs on a specific rank
        weight_column='corrected_weights',   # options: weights, corrected_weights, bitscore, score, weighted_score, corrected_weighted_score
        protein_column="sseqid",             # name of column containing protein accessions, will be retained according to lca 
        taxid_column="OX",                   # name of column containing taxonomy ids, will be retained according to lca 
        weight_cutoff=0.6,                   # minimum fraction of total weights 
        
        filter_cutoff=5,                     # post-lca filtering                  
         
        write_database="Proteins",           # Options (False, "Proteins","Taxids"): do not make a database(False), use aligned proteins ("Proteins") use aligned taxids("Taxids").
        proteins=""
        
        ):

    # #%%
    # #test
    # df=target_decoys[0]                                  #dataframe with at least a column called Peptide, and ranks
    
    # Output_directory
    # group_on="Peptide"                  
    # denovo_peptides=None                #original PSM DataFrame, used to merge back results
    
    # #Weighing parameters
    # method="weighted"                   # Options: False (conventional LCA), "weighted" (weighted lca), "focused" (focused lca)
    # weight_rank="genus"                 # weighted weighs on a specific rank
    # weight_column='corrected_weights'   # options: weights, corrected_weights, bitscore, score, weighted_score, corrected_weighted_score
    # protein_column="sseqid"            # name of column containing protein accessions, will be retained according to lca 
    # taxid_column="OX"
    # weight_cutoff=0.6                   # minimum fraction of total weights 
    
    # filter_cutoff=5                     # post-lca filtering                  
     
    # write_database=True        
    # proteins=""
    
    if method=="focused":
        lin=[]
        for rank in ranks:
            wc=df.groupby([group_on,rank]).agg({weight_column:['sum']})/df.groupby(group_on).agg({weight_column:['sum']})
            wc.columns=wc.columns.droplevel(1)
            wc=wc.reset_index(rank)
            lin.append(wc[wc[weight_column]>=weight_cutoff][rank])
            
        lin=pd.concat(lin,axis=1)
        lcas=pd.DataFrame(df[group_on]).drop_duplicates().merge(lin,on=group_on).set_index(group_on)
        
    else:
    
        if method=="weighted":
            wc=df.groupby([group_on,weight_rank]).agg({weight_column:['sum']})/df.groupby(group_on).agg({weight_column:['sum']})
            wc.columns=wc.columns.droplevel(1)
            group=wc.reset_index().sort_values(by=[group_on,weight_column],ascending=False).groupby(group_on,sort=False)
            cut_df=pd.concat([x.iloc[0:(x[weight_column].cumsum()<weight_cutoff).sum()+1] for n,x in group])
            df=df[(df[group_on]+df[weight_rank]).isin(cut_df[group_on]+cut_df[weight_rank])]
        
        lcas=df.groupby(group_on)[ranks].nth(0)[(df.groupby(group_on)[ranks].nunique()==1)] #vectorized standard lca
        
    #aggregate proteins and taxids that follow lca:
    last=lcas.fillna(method="ffill",axis=1).iloc[:,-1]
    if protein_column in df.columns:
        lcas["proteins"]=df[df[ranks].add(df[group_on],axis=0).isin(last.tolist()+last.index).any(axis=1)].groupby(group_on)[protein_column].apply(lambda x: ", ".join(list(set(x))))
    lcas["taxids"]=df[df[ranks].add(df[group_on],axis=0).isin(last.tolist()+last.index).any(axis=1)].groupby(group_on)[taxid_column].apply(lambda x: ", ".join(list(set(x))))

    #add back proteins with no common ancestor concensus
    no_lca=df[~df[group_on].isin(lcas.index)]
    no_lca=pd.DataFrame(no_lca.groupby(group_on)[protein_column].apply(lambda x: ", ".join(x)))
    no_lca.columns=["proteins"]
    no_lca[ranks]=[""]*len(ranks)
    lcas=pd.concat([lcas,no_lca],axis=0).fillna("")

    #write outputs (lca and de novo peptides)
    name_dict={"weights":"w","corrected_weights":"cw","bitscore":"b","score":"s","weighted_score":"ws","corrected_weighted_score":"cws"}
    outdirs=["lca"]
    
    if type(denovo_peptides)!=type(None): outdirs+=["PSMs"]
    if write_database:                    outdirs+=["Database"]
    [os.mkdir(str(Path(Output_directory,d))) for d in outdirs if not os.path.exists(str(Path(Output_directory,d)))]

    prefix=""
    if method=="weighted": prefix="weighted_"+name_dict.get(weight_column)+str(weight_cutoff).replace(".","_")+"_" 
    if method=="focused":  prefix="focused_" +name_dict.get(weight_column)+str(weight_cutoff).replace(".","_")+"_" 
    
    #postfilter
    removed_taxons=[]
    if type(denovo_peptides)!=type(None):
        denovo_peptides,lcas,removed_taxons=Postfilter(lcas,denovo_peptides=denovo_peptides,group_on=group_on,filter_cutoff=filter_cutoff) 
        
        PSMs_out=str(Path(Output_directory,"PSMs",prefix+"PSMs.tsv"))
        denovo_peptides.to_csv(PSMs_out,sep="\t") #write_PSMs
            
    #filter proteins and taxids with postfilter 
    fdf=df[df[ranks].isin(removed_taxons).any(axis=1)]
    if protein_column in fdf.columns:
        removed_proteins=fdf[protein_column]
        exProts=lcas.explode("proteins")["proteins"]
        lcas["proteins"]=exProts[~exProts["proteins"].isin(removed_proteins).groupby(exProts.index).apply(list) 
    exOX=lcas.explode("taxids")
    allowed_taxids=exOX[~exOX["taxids"].isin(removed_taxons)
    lcas["taxids"]=allowed_taxids.groupby(exOX.index).apply(list)
                                  
    lcas.to_csv(str(Path(Output_directory,"lca",prefix+"lca.tsv")),sep="\t") #write_lca
    
    #write_database 
    db_out=str(Path(Output_directory,"Database",prefix+"Database.fa"))
    if write_database=="Proteins":
        fasta_str=Proteins_to_fasta(df,lcas,removed_proteins=removed_proteins)
        if fasta_str:
            with open(db_out,"w") as f:
                f.write(fasta_str)
    if write_database=="Taxids":
        Taxids_to_fasta(lineages=allowed_taxids.tolist())
              
                
    if type(denovo_peptides)!=type(None):
            
        return PSMs_out


