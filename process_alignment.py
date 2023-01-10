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
from __main__ import ranks, ncbi_taxdf, de_novo_file
import re

# read alignment file to generator and filter on topx % scoring
def Diamond_alignment_Reader(input_file,BIT,Output_directory,
                             filter_dynamic_score=True,
                             filter_alc_bit=False,
                             score_cutoff=0.9,
                             Read_batch=1000000,
                             output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]):
    
    cdf=pd.read_csv(input_file, sep='\t', chunksize=Read_batch,names=output_columns) # read to generator
 
    sc=[]
    for ix,c in enumerate(cdf):
        _,index=np.unique(c.qseqid,return_index=True)
        index.sort()
        index_m=index.tolist()+[max(index)+1]
       
        if ix>0:
            c=pd.concat([sc[-1],c])
       
        sc=[c.iloc[index_m[n]:index_m[n+1]] for n in range(len(index_m)-1)]
        for s in sc[:-1]:
            s=s.assign(Peptide=s.qseqid.apply(lambda x: x.split(";")[0]))
            if filter_dynamic_score:
                s.loc[:,"qlen"]=s["Peptide"].apply(len)
                s.loc[:,"slen"]=s["sseq"].apply(len)
                s.loc[:,"score"]=s["bitscore"]*s["slen"]/s["qlen"]*s["pident"]/100
                s=s[s.score>=s.score.max()*(score_cutoff)]
            elif filter_alc_bit:
                s.loc[:,"qlen"]=s["Peptide"].apply(len)
                s.loc[:,"slen"]=s["sseq"].apply(len)
                s.loc[:,"score"]=s["bitscore"]*s["slen"]/s["qlen"]*s["pident"]/100
                s=s[s.bitscore>=BIT]
            yield s
    
    # last one
    s=sc[-1]
    s=s.assign(Peptide=s.qseqid.apply(lambda x: x.split(";")[0]))
    if filter_dynamic_score:
        s.loc[:,"qlen"]=s["Peptide"].apply(len)
        s.loc[:,"slen"]=s["sseq"].apply(len)
        s.loc[:,"score"]=s["bitscore"]*s["slen"]/s["qlen"]*s["pident"]/100
        s=s[s.score>=s.score.max()*(score_cutoff)] 
    elif filter_alc_bit:
        s.loc[:,"qlen"]=s["Peptide"].apply(len)
        s.loc[:,"slen"]=s["sseq"].apply(len)
        s.loc[:,"score"]=s["bitscore"]*s["slen"]/s["qlen"]*s["pident"]/100
        s=s[s.bitscore>=BIT]
    yield s
 

def fill_g(x):
    cons=[np.array(list(g)) for g in mit.consecutive_groups(np.argwhere(x==""))]
    if cons:
        for c in cons:  
            c=c.flatten().tolist()
            if c[-1]!=len(x)-1:
                for i in c:
                    x[i]="gap_"+str(x[c[-1]+1])+"_"+ranks[i]
    return list(x)


def Process_diamond_alignment(alignments,BIT,Output_directory,
                              filter_dynamic_score=True,
                              filter_alc_bit=False,
                              fill_gaps=True,
                              Target_Decoy_prefilter=False,
                              prefilter_quantile=0.95,
                              weight_rank="genus"):
    
    target_decoys=[]
    for alignment in alignments:
        
        iterable=Diamond_alignment_Reader(alignment,BIT,Output_directory,filter_dynamic_score,filter_alc_bit)
        al=pd.concat([i for i in iterable])
       
        # add NCBI taxonomy UniRef or UniprotKB headers
        if  re.search("UniRef", al['sseqid'].iloc[1]):
            taxids=pd.DataFrame(al["stitle"].apply(lambda x: ((x+"TaxID=").split("TaxID=")[1].split()+[" "])[0]))
            taxids=taxids.replace('N/A','131567') # replace missing taxa with "cellular organism"
            taxids['stitle']=pd.to_numeric(taxids['stitle'])
            al["OX"]=taxids     
            al["gene"]=al["stitle"].apply(lambda x: ((x+"RepID=").split("RepID=")[1].split()+[" "])[0])
        
        else:
            al["OX"]=al["stitle"].apply(lambda x: ((x+"OX= ").split("OX=")[1].split()+[" "])[0]).astype(int)
            al["gene"]=al["stitle"].apply(lambda x: ((x+"GN= ").split("GN=")[1].split()+[" "])[0])
        
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
                mw["corrected_weights"]=((mw.ntarget_count-mw.ndecoy_count)/(mw.ntarget_count+mw.ndecoy_count)+1)/2           # percentile difference normalized to 0-1
                mw["ntarget_count"]=(mw["ntarget_count"]-mw["ntarget_count"].min())/(mw["ntarget_count"].max()-mw["ntarget_count"].min()) #renormalize to 0-1
                mw["ndecoy_count"]=(mw["ndecoy_count"]-mw["ndecoy_count"].min())/(mw["ndecoy_count"].max()-mw["ndecoy_count"].min()) #renormalize to 0-1
                target=target.merge(mw[[weight_rank,"ntarget_count","corrected_weights"]].rename(columns={"ntarget_count":"weights"}),on=weight_rank)
                decoy=  decoy.merge(mw[[weight_rank, "ndecoy_count","corrected_weights"]].rename(columns={ "ndecoy_count":"weights"}),on=weight_rank)

                target["weighted_score"]=target["score"]*target["weights"]
                target["corrected_weighted_score"]=target["score"]*target["corrected_weights"]
                decoy["weighted_score"]=decoy["score"]* decoy["weights"]
                decoy["corrected_weighted_score"]=decoy["score"]* decoy["corrected_weights"]
            
            if alignment.find('only')!=-1:
                target=target.replace('alignment','DN_only_alignment', regex=True)
                decoy=decoy.replace('alignment','DN_only_alignment', regex=True)
                
            target_decoys.extend([target,decoy])
            
        else:
            target=al
            decoy=""
            target_decoys.extend([target,decoy])
 
    return target_decoys


def lca(df,
        denovo_peptides="", 
        Output_directory="", 
        group_on="Peptide",
        weighted=False,
        weight_rank="genus",
        weight_column="weights",   # "weights" for weighted lca, "corrected_weights" for corrected_weighted lca
        weight_cutoff=0.6   # minimum fraction of total weights
        ):

    # define filename prefix
    if weighted==True:
        if weight_column=="weights":
            T="_w"
        else:
            T="_cw"
    else:
        T="_"

    # start lca
    x=df.drop_duplicates().set_index(group_on).groupby(group_on)             
    lcas=[]
    for ix,g in enumerate(x):  
        q=  [g[0]]
        group=g[1]
        lca=[]

        if len(group)>1:
           
            if weighted==True:
                # allow only taxa with a weight > cut off
                if len(group[group[weight_rank]!=""])>1:
                    group=group[group[weight_rank]!=""]
                    w=group.groupby(weight_rank)[weight_column].sum().sort_values(ascending=False)
                    group=group[group[weight_rank].isin(w.index[0:(w.cumsum()/w.sum()<weight_cutoff).sum()+1])]    
    
            for ir,rank in enumerate(ranks):
                if group[rank].nunique()!=1:
                    break
   
            if ir==len(ranks)+1:
                ir-=1
    
            lca=group.iloc[0][ranks[0:ir]].tolist()
            lca+=[""]*(len(ranks)-len(lca))
     
            lcas.append(q+lca)
     
        elif len(group)==1:
            lca=group.iloc[0][ranks].tolist()
            lcas.append(q+lca)

    ldf=pd.DataFrame(lcas,columns=[g[1].index.name]+ranks)
    
    if "Peptide" in ldf.columns:
        
        # check if output folder exists
        output_folder=str(Path(Output_directory,"lca"))
        isExist = os.path.exists(output_folder)
        if not isExist:
              os.mkdir(output_folder)

        # check if target/decoy alignment file
        if  re.search(":'Target'", df['qseqid'].iloc[1]):
            
            if re.search("DN_only", df['qseqid'].iloc[1]):
                # add taxonomies to denovo file
                denovo_peptides_lca=denovo_peptides.merge(ldf,on="Peptide",how="left")
                denovo_peptides_lca=Postfilter(denovo_peptides_lca)
                # write results
                denovo_peptides_lca.to_csv(str(Path(output_folder,Path(de_novo_file).stem+T+"lca_DN_only_combined.tsv")),sep="\t")
            else:
                # add taxonomies to denovo file
                denovo_peptides_lca=denovo_peptides.merge(ldf,on="Peptide",how="left")
                denovo_peptides_lca=Postfilter(denovo_peptides_lca)
                # write results
                denovo_peptides_lca.to_csv(str(Path(output_folder,Path(de_novo_file).stem+T+"lca_combined.tsv")),sep="\t")
                return denovo_peptides_lca

        else:
            # write results
            ldf=Postfilter(ldf)
            ldf.to_csv(str(Path(output_folder,Path(de_novo_file).stem+T+"lca_decoy_combined.tsv")),sep="\t")

    else:
        return ldf # only for DB searching outputs


def score_lca(df,
              denovo_peptides="", 
              Output_directory="", 
              score_column="weighted_score",    # bitscore, weighed_score, corrected_weighted_score
              group_on="Peptide",
              score_cutoff=0.6,         # minimum fraction of total score,
              rescoring=False,          # if true, total bitscore, recalculated every taxonomic rank
              ): 
    
              # start lca
              x=df.drop_duplicates().set_index(group_on).groupby(group_on) 
              lcas=[]       

              for ix,g in enumerate(x):
                  q=  [g[0]]
                  group=g[1]
                  t=[]
                 
                  sbit=group[score_column].sum()*score_cutoff # minimum fraction of score that a taxonomy should have
                  if len(group)>1:
                      lca=[]
                      for ir,rank in enumerate(ranks):
                                                       
                          s=group.groupby(rank)[score_column].sum().sort_values()
            
                          if s.max()>sbit:
                              t=s.index[-1]
                              lca.append(t)
                              
                              if rescoring==True:
                                  group=group[group[rank]==t] # keep only dominant taxon
                                  sbit=group[score_column].sum()*score_cutoff # minimum fraction of score that a taxonomy should have
                          else:
                              ir-=1
                              break
                       
                          if ir==len(ranks)+1:
                              ir-=1
                              
                      tax=lca+[""]*(len(ranks)-len(lca))
                      lcas.append(q+tax)
          
                  elif len(group)==1:
                      lca=group.iloc[0][ranks].tolist()
                      lcas.append(q+lca)
 
              ldf=pd.DataFrame(lcas,columns=[g[1].index.name]+ranks)

              if "Peptide" in ldf.columns:
                  
                  # check if output folder exists
                  output_folder=str(Path(Output_directory,"lca"))
                  isExist = os.path.exists(output_folder)
                  if not isExist:
                     os.mkdir(output_folder)
                
                  # check if target/decoy alignment file
                  if  re.search(":'Target'", df['qseqid'].iloc[1]):
                      
                      if re.search("DN_only", df['qseqid'].iloc[1]):
                          # add taxonomies to denovo file
                          denovo_peptides_blca=denovo_peptides.merge(ldf,on="Peptide",how="left")
                          denovo_peptides_blca=Postfilter(denovo_peptides_blca)
                          # write results
                          denovo_peptides_blca.to_csv(str(Path(output_folder,Path(de_novo_file).stem+"_blca_DN_only_combined.tsv")),sep="\t")
                      else:
                          # add taxonomies to denovo file
                          denovo_peptides_blca=denovo_peptides.merge(ldf,on="Peptide",how="left")
                          denovo_peptides_blca=Postfilter(denovo_peptides_blca)
                          # write results
                          denovo_peptides_blca.to_csv(str(Path(output_folder,Path(de_novo_file).stem+"_blca_combined.tsv")),sep="\t")
                          return denovo_peptides_blca
                  
                  else:
                     # write results
                     ldf=Postfilter(ldf)
                     ldf.to_csv(str(Path(output_folder,Path(de_novo_file).stem+"_blca_decoy_combined.tsv")),sep="\t")

              else:
                return ldf

# Postfilter
def Postfilter(tpeptide_df, # target peptides
               dpeptide_df=None, # decoy peptides 
               weight_rank="genus",
               fixed_cutoff=5, # fixed branch cutoff / as used in tax grouping code
                ):
 
    # fixed denoising
    if dpeptide_df==None:
        
        unirows=tpeptide_df[["Peptide"]+ranks].drop_duplicates()[ranks].astype(str).values.tolist()
        jbranch=["#*".join(i) for i in unirows]
        fbranch=[branch for branch, counts in Counter(jbranch).items() if counts >= fixed_cutoff]
        allowed_taxids=set(chain(*[i.split("#*") for i in fbranch]))       
        for i in ranks:
            tpeptide_df.loc[~tpeptide_df[i].isin(allowed_taxids),i]=""
            
        return tpeptide_df
    
