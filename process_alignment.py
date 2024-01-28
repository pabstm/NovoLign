# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 15:17:18 2024

@author: e_kle
"""

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
from __main__ import ranks, ncbi_taxdf, fasta_database_path
from Bio import SeqIO
import itertools

import warnings
warnings.filterwarnings("ignore")

# read alignment file and filter for bit score and ALC
def filter_chunk(s,top_score_fraction,minimum_bitscore,min_ALC_score):
    
    s=s.assign(Peptide=s.qseqid.apply(lambda x: x.split(";")[0]))
    s.loc[:,"qlen"]=s["Peptide"].apply(len)
    s.loc[:,"slen"]=s["sseq"].apply(len)
    s.loc[:,"score"]=s["bitscore"]*s["slen"]/s["qlen"]*s["pident"]/100 #alternative score which considers alignment length
    
    #filter based on bitscore
    if minimum_bitscore:      s=s[s.bitscore>=minimum_bitscore]
    
    if s.empty == False:
        #extract and filter based on ALC
        split_list=[]
        for index, row in s.iterrows():
            split_list.append(row["qseqid"].split("'ALC(%)':")[1])
        split_list=[float(item.replace("'","").replace("}","")) 
                      for item in split_list]
        alc=pd.DataFrame(split_list).set_index(s.index)
        alc.rename(columns={0:"ALC"}, inplace=True)
        if min_ALC_score:         s=s[alc.ALC>=min_ALC_score]
    
    return s


def Diamond_alignment_Reader(input_file,
                             Output_directory,
                             top_score_fraction=0.9,
                             minimum_bitscore=0,
                             min_ALC_score=40,
                             Read_batch=1000000,
                             output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]):

    cdf=pd.read_csv(input_file, sep='\t', chunksize=Read_batch,names=output_columns)
 
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
            yield filter_chunk(s,top_score_fraction=top_score_fraction,minimum_bitscore=minimum_bitscore,min_ALC_score=min_ALC_score)
        
    s=sc[-1]
    yield filter_chunk(s,top_score_fraction=top_score_fraction,minimum_bitscore=minimum_bitscore,min_ALC_score=min_ALC_score)
 

def fill_g(x):
    cons=[np.array(list(g)) for g in mit.consecutive_groups(np.argwhere(x==""))]
    if cons:
        for c in cons:  
            c=c.flatten().tolist()
            if c[-1]!=len(x)-1:
                for i in c:
                    x[i]="gap_"+str(x[c[-1]+1])+"_"+ranks[i]
    return list(x)


def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g
        

def Process_diamond_alignment(alignments,
                              Output_directory,
                              fill_gaps=True,
                              #fill_gaps=False,
                              top_score_fraction=0.9,
                              minimum_bitscore=0,
                              min_ALC_score=40,
                              weight_rank="genus",
                              output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]
                              ):
   
    target_decoys=[]
    if type(alignments)==str: alignments=[alignments]
    for alignment in alignments:
       
        iterable=Diamond_alignment_Reader(alignment,Output_directory,top_score_fraction=top_score_fraction,minimum_bitscore=minimum_bitscore,min_ALC_score=min_ALC_score,output_columns=output_columns)
        al=pd.concat([i for i in iterable])
       
        #add NCBI taxonomy for UniRef or UniprotKB headers
        if al.head(10)["stitle"].str.contains("TaxID=").any(): 
            al["OX"] = al["stitle"].apply(lambda x: ((x+"TaxID="+" ").split("TaxID=")[1].split()+[" "])[0].replace("N/A","131567")).astype(int)
            al["gene"]=al["stitle"].apply(lambda x: ((x+"RepID="+" ").split("RepID=")[1].split()+[" "])[0])
        else:
            al["OX"]  =al["stitle"].apply(lambda x: ((x+"OX="+" ").split("OX=")[1].split()+[" "])[0].replace("N/A","131567")).astype(int)
            al["gene"]=al["stitle"].apply(lambda x: ((x+"GN="+" ").split("GN=")[1].split()+[" "])[0])

        al=al.merge(ncbi_taxdf,how="left",on="OX").fillna("")
   

        #fill annotation gaps in lineages
        if fill_gaps:
            gaps=al[(al[ranks]=="").any(axis=1)]
            u,ix,inv=np.unique(gaps[ranks].apply(";".join,axis=1),return_inverse=True,return_index=True)
            u=gaps.iloc[ix][ranks].values
            gap_1=pd.DataFrame(np.array(list((map(fill_g,u))))[inv]).set_index(gaps.index)
            al.loc[gap_1.index,ranks]=gap_1.values

        al=pd.concat([al, #add information stored in the header of the fasta file to the read dataframe
                  pd.DataFrame.from_dict(al.qseqid.str.split(";").apply(lambda x: eval(x[1])).values.tolist(),orient="columns")],axis=1)



        if "alignment_Target_Decoy" in al.columns:
            target=al[al["alignment_Target_Decoy"]=="Target"]
            decoy=al[al["alignment_Target_Decoy"]=="Decoy"]

            #add weights
            if type(weight_rank)!=type(None):

                #add weights
                mw=pd.concat([d.groupby(weight_rank).size() for d in [target,decoy]],axis=1).fillna(0).reset_index()
                mw.columns=[weight_rank,"target_count","decoy_count"]
                mw["ntarget_count"]=mw["target_count"]/mw["target_count"].sum()
                mw["ndecoy_count"]=mw["decoy_count"]/mw["decoy_count"].sum()
                mw["corrected_weights"]=((mw.ntarget_count-mw.ndecoy_count)/(mw.ntarget_count+mw.ndecoy_count)+1)/2                        # percentile difference normalized to 0-1
                mw["ntarget_count"]=(mw["ntarget_count"]-mw["ntarget_count"].min())/(mw["ntarget_count"].max()-mw["ntarget_count"].min())  # renormalize to 0-1
                mw["ndecoy_count"]=(mw["ndecoy_count"]-mw["ndecoy_count"].min())/(mw["ndecoy_count"].max()-mw["ndecoy_count"].min())       # renormalize to 0-1
                target=target.merge(mw[[weight_rank,"ntarget_count","corrected_weights"]].rename(columns={"ntarget_count":"weights"}),on=weight_rank)
                decoy=decoy.merge(mw[[weight_rank, "ndecoy_count","corrected_weights"]].rename(columns={ "ndecoy_count":"weights"}),on=weight_rank)
                target["weighted_score"]=target["score"]*target["weights"]
                target["corrected_weighted_score"]=target["score"]*target["corrected_weights"]
                decoy["weighted_score"]=decoy["score"]* decoy["weights"]
                decoy["corrected_weighted_score"]=decoy["score"]* decoy["corrected_weights"]
                
                target_decoys.append(pd.concat([target,decoy]))
         
        else:
            target_decoys.append(al)

    return target_decoys

def Postfilter(lcas,
               denovo_peptides=None,
               group_on="Scan",      
               filter_cutoff=5,                    
                                ): 
  
    pdf=lcas
    unirows=pdf[ranks].groupby(ranks).size().rename("Count").reset_index()           
    allowed_taxons=np.unique(unirows.loc[unirows["Count"]>filter_cutoff,ranks])
    removed_taxons=np.unique(unirows.loc[unirows["Count"]<=filter_cutoff,ranks])
    removed_taxons=set(removed_taxons)-set(allowed_taxons)
    
    for i in ranks:
        lcas.loc[~lcas[i].isin(allowed_taxons),i]=""

    return lcas,removed_taxons
    

def Taxids_to_fasta(lineages,                       #list of taxons or taxonomy ids
                    output_path,                    #File to which the resulting database is written
                    Database=fasta_database_path,   #Database from which to select the sequences from 
                          #
                    minimum_rank="OX",              #Minimum rank that supplied taxonomic lineages should have
                    add_decoy=False,                #Add reverse decoy of aligned proteins to database
                    alignment_df=[],                #Used to add the de novo peptides to the database
                        ):
    
        #This code takes in a list of lineages and checks whether the lineages in it appear in the NCBI Taxonomy Database. 
        #If they do, it will get the unique taxonomic ranks for the lineages starting at a given minimum rank.
        #The set of lineages is then expanded with the unique Taxonomic Identification numbers and then stored in a list.
        lineages=[str(i) for i in lineages if i !=""]
        lineages=[int(i)  if i.isdigit() else i for i in list(set(lineages))]
        rank_order=np.array(ranks+["OS","OX"])
        allowed_ranks=rank_order[np.argwhere(rank_order==minimum_rank)[0][0]:]
        expanded_lineages=ncbi_taxdf.loc[ncbi_taxdf[allowed_ranks].isin(lineages).any(axis=1),allowed_ranks]
        expanded_taxids=list(set(ncbi_taxdf.loc[ncbi_taxdf.isin(expanded_lineages.values.flatten()).any(axis=1),"OX"]))
  
        #determine the taxonomic rank of each input
        recs=SeqIO.parse(Database,format="fasta")
        chunks=chunk_gen(recs)    
        
        with open(output_path,"w") as f:
            pass
        with open(output_path,"a") as f:
            
            for ic,c in enumerate(chunks):
                print("chunk "+ str(ic))
        
                chunk_df=pd.DataFrame([[str(r.seq),r.description,r.id] for r in c],columns=["seq","description","id"])
                taxids=chunk_df["description"].str.split("OX=").apply(lambda x: x[-1].split(" ")[0]).astype(int)
                chunk_df=chunk_df[taxids.isin(expanded_taxids)]
                f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")
                if add_decoy:
                    f.write("\n".join(">decoy_"+chunk_df["description"]+"\n"+chunk_df["seq"].str[::-1])+"\n")
                                                                                                 
            if len(alignment_df):
                dn_peptides=alignment_df.qseq.drop_duplicates().reset_index()
                f.write("\n".join((">"+"DeNovo_"+dn_peptides["index"].astype(str)+"\n"+dn_peptides["qseq"]))+"\n")                                                                                                 
    
        return output_path

#%%   
def Proteins_to_fasta(df,                       #processed diamond alignment 
                      lcas=None,                #input: lca df, used to generate list of proteins to include in database
                      removed_proteins=[],      #input: list of proteins that are excluded from database
                      minimum_rank="",          #input: taxonomic rank, in case lca is supplied, this can only select lcas with a certain level of taxonomic detail 
                      add_decoy=False,          #add reverse decoy of aligned proteins to database
                      add_denovo_peptides=True, #add de novo peptides to database (not included in decoy)
                      ):
  
    if "full_sseq" not in df.columns:
        print("No full protein sequences found, run diamond alignment with full_sseq in output columns!")
        return
    
    #filter database based on lca     
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
  
    return fasta_str


def lca(df,                                  # dataframe with columns named Peptide/Scan and ranks
        
        Output_directory,
        group_on="Scan",                  
        denovo_peptides=None,                # crude PSM dataframe

        #Weighing parameters
        method="weighted",                   # Options: False (conventional LCA), "weighted" (weighted lca), "bitscore" (bitscore lca)
        weight_rank="genus",                 # weighted weighs on a specific rank
        weight_column='corrected_weights',   # options: weights, corrected_weights, bitscore, score, weighted_score, corrected_weighted_score
        protein_column="sseqid",             # name of column containing protein accessions, will be retained according to lca 
        taxid_column="OX",                   # name of column containing taxonomy ids, will be retained according to lca 
        weight_cutoff=0.6,                   # minimum fraction of total weights 
        filter_cutoff=5,                     # post-lca filtering                  
        
        #options for database construction
        write_database="False",              # Options ("False","Proteins","Taxids"): do not make a database(False), use aligned proteins ("Proteins") use aligned taxids("Taxids").
        minimum_rank="OX",                   # Only used when writing database based on taxonomy
        add_decoy=False,                     # add a reversed decoy sequences to the database
        add_denovo=False,                    # add the denovo sequences peptides to the database
        ):
    
    #separate alignments based on target/decoys
    if df.columns.str.contains('alignment_Target_Decoy').any():
        df_target = df.loc[df['alignment_Target_Decoy'] == 'Target']
        df_decoy = df.loc[df['alignment_Target_Decoy'] == 'Decoy']
        if len(df_decoy.index)==0:
            df_decoy.loc[0] = df_target.loc[0]
            df_decoy.loc[0,'alignment_Target_Decoy'] = 'Decoy' 
        df_list = [df_decoy,df_target]
    else:
        df_target = df
        df_list = [df_target]
    
    #loop through target/decoys
    for df in df_list:
    
        if Output_directory!="simple_lca":
            #df[taxid_column]=df[taxid_column].astype(str)
            df.loc[:,taxid_column] = df[taxid_column].astype(str)
        
        if method=="bitscore":
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
            
            lcas=df.groupby(group_on)[ranks].nth(0)[(df.groupby(group_on)[ranks].nunique()==1)] #vectorized conventional lca
            
            if Output_directory == "simple_lca":
                return lcas
                
        #aggregate proteins and taxids that follow lca
        last=lcas.fillna(method="ffill",axis=1).iloc[:,-1]
        if protein_column in df.columns:
            lcas["proteins"]=df[df[ranks].add(df[group_on],axis=0).isin(last.tolist()+last.index).any(axis=1)].groupby(group_on)[protein_column].apply(lambda x: list(set(x)))
        lcas["taxids"]=df[df[ranks].add(df[group_on],axis=0).isin(last.tolist()+last.index).any(axis=1)].groupby(group_on)[taxid_column].apply(lambda x: list(set(x)))
        
        #add back proteins without common ancestor concensus
        no_lca=df[~df[group_on].isin(lcas.index)]
        no_lca=pd.DataFrame(no_lca.groupby(group_on)[protein_column].apply(lambda x: ", ".join(x)))
        no_lca.columns=["proteins"]
        no_lca[ranks]=[""]*len(ranks)
        lcas=pd.concat([lcas,no_lca],axis=0).fillna("")
    
        #write outputs
        name_dict={"weights":"W","corrected_weights":"CW","bitscore":"BIT","score":"S","weighted_score":"WS","corrected_weighted_score":"CWS"}
        outdirs=["lca"]
        
        if type(denovo_peptides)!=type(None): outdirs+=["psms"]
        if write_database:                    outdirs+=["database"]
        [os.mkdir(str(Path(Output_directory,d))) for d in outdirs if not os.path.exists(str(Path(Output_directory,d)))]
    
        prefix=""
        if method=="weighted": prefix=name_dict.get(weight_column)+"_"+str(weight_cutoff).replace(".","")+"_" 
        if method=="bitscore":  prefix=name_dict.get(weight_column)+"_"+str(weight_cutoff).replace(".","")+"_" 
        if method=="standard": prefix="CON_0_"

        #filter proteins and taxids with postfilter 
        if filter_cutoff:
            lcas,removed_taxons=Postfilter(lcas,denovo_peptides=denovo_peptides,group_on=group_on,filter_cutoff=filter_cutoff) 
            
            fdf=df[df[ranks].isin(removed_taxons).any(axis=1)]
            if protein_column in fdf.columns:
                removed_proteins=fdf[protein_column]
                exProts=lcas.explode("proteins")["proteins"]
                exProts=exProts[~exProts.isin(removed_proteins)]
                lcas["proteins"]=exProts.groupby(exProts.index).apply(list)
                
            removed_taxids=fdf[taxid_column]
            exOX=lcas.explode("taxids")["taxids"]
            exOX=exOX[~exOX.isin(removed_taxids)]
            lcas["taxids"]=exOX.groupby(exOX.index).apply(list)
    
        allowed_taxids=np.unique(lcas.explode("taxids")["taxids"].dropna())
    
        lcas.loc[lcas["proteins"].isnull(),"proteins"] = lcas.loc[lcas["proteins"].isnull(),"proteins"].apply(lambda x: [""])
        lcas.loc[lcas["taxids"].isnull(),"taxids"] = lcas.loc[lcas["taxids"].isnull(),"taxids"].apply(lambda x: [""])
        lcas["proteins"]=lcas["proteins"].apply(lambda x: ", ".join(x))  
        lcas["taxids"]=lcas["taxids"].apply(lambda x: ", ".join(x))
        
        if df.iloc[0, df.columns.get_loc('alignment_Target_Decoy')] == 'Decoy':
            file=str(Path(Output_directory,"lca",prefix+str(filter_cutoff)+"_DECOY_lca.tsv"))
            lcas.to_csv(file,sep="\t")
        else:
            file=str(Path(Output_directory,"lca",prefix+str(filter_cutoff)+"_TARGET_lca.tsv")) 
            lcas.to_csv(file,sep="\t")            
          
            #write DN bitscore database
            db_out=str(Path(Output_directory,"database",prefix+write_database+"_"+minimum_rank+"_DN_bitscore_referenc_database.fasta"))
            if write_database=="Proteins":
                fasta_str=Proteins_to_fasta(df,lcas,removed_proteins=removed_proteins,add_denovo=add_denovo)
                if fasta_str:
                    with open(db_out,"w") as f:
                        f.write(fasta_str)
            if write_database=="Taxids":
                if add_denovo:
                    Taxids_to_fasta(lineages=allowed_taxids,Database=fasta_database_path,alignment_df=df,output_path=db_out,minimum_rank=minimum_rank,add_decoy=add_decoy)
                else:
                    Taxids_to_fasta(lineages=allowed_taxids,Database=fasta_database_path,output_path=db_out,minimum_rank=minimum_rank,add_decoy=add_decoy)
                      
            if type(denovo_peptides)!=type(None):
                PSMs_out=str(Path(Output_directory,"psms",prefix+"PSMs.tsv"))
                denovo_peptides.to_csv(PSMs_out,sep="\t")

            return file
        
        
