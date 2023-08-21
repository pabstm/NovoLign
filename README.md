# NovoLign



This is the repository for the NovoLign pipeline, as described in:<br>
"NovoLign: metaproteomic profiling by de novo sequence alignment" 


The pipeline was established and tested with shotgun (meta)proteomics data obtained from Q Exactive Orbitrap Mass Spectrometers, using PEAKS generated de novo sequence lists. The generation of accurate de novo peptide sequence lists depends on high quality peptide sequencing spectra. NovoLign has been tested and developed in an Anaconda Spyder environment.

<img src="https://github.com/hbckleikamp/NovoLign/blob/main/images/workflow.svg" width="300" height="400" align="right">
<br>

#### What is NovoLign? 
NovoLign is an tool that uses DIAMOND for high-througput, error-tolerant annotation of de novo sequenced peptides.
As de novo sequencing is independent of database composition, it provides an unbiased alternative to conventional database searching methods. By aligning de novo peptides to general databases such as UniRef100, NovoLign can find related peptide sequences, which can aid in construction of sample specific databases, database and experiment quality control, and asses the completeness of references databases.


## Basic use

#### Setting up NovoLign 
NovoLign is designed to interface with multiple commonly used reference databases, belonging to NCBI, UniprotKB or GTDB.
Setup scripts are supplied in the folder `Setup` and can be run automatically with `run_Setup.py`
To setup protein databases for GTDB, utility scripts are supplied to download unzip and merge GTDB reference proteomes, and add taxonomies to fasta headers. For NCBI, a linearized taxonomy is constructed for specified taxonomic ranks using taxdump.
scripts are supplied for database curation, such as euqating I and L, and removing ambiguous amino acids.
Lastly, scripts are supplied for automated downloading of DIAMOND.

<br>

#### How does it work? (Placeholder)

The NovoLign pipeline consists of 5 parts:
1. *DIAMOND alignment* 
2. *Lowest Common Ancestor analysis (LCA)*
3. *Taxonomy report*
4. *Spectral quality report*
5. *Database coverage report*

In Part 1: input files are read, parsed, filtered and submitted to Unipept for taxonomic and functional annotations. <br>
In Part 2: Unipept taxonomic annotations are quantified, and visual outputs are generated. <br>
In Part 3: Unipept functional annotations are matched to KEGG orthologies and quantified. <br>
In Part 4: Unipept functional annotations are matched to KEGG orthologies and quantified. <br>
In Part 5: Unipept functional annotations are matched to KEGG orthologies and quantified. <br>

<br>
Using alignment parameters optimized to short peptide homology

LCAs optimized for short peptide homology

#### What does it do? (Placeholder)
<br>
Outputs


   # 1. Conventional lca
        denovo_peptides_lca=lca(target_decoy,Output_directory,denovo_peptides=denovo_peptides,method="standard",filter_cutoff=freq_cut,minimum_rank=DB_rank,write_database=DB)
        # 2. Bitsore lca 
        denovo_peptides_blca=lca(target_decoy,Output_directory,denovo_peptides=denovo_peptides,method="focused",weight_column="bitscore",filter_cutoff=freq_cut)
        # # 3. Weighted lca
        denovo_peptides_wlca=lca(target_decoy,Output_directory,denovo_peptides=denovo_peptides,method="weighted",weight_column="weights",filter_cutoff=freq_cut)
     

#### Running NovoLign 
- Novobridge is designed as a single "tunable" python script.
- Novobridge does not offer command line options, but parameters can be altered in the script Novobridge.py



#### What input files does it use? 

<br>

#### What outputs does it generate? 

<br>

## Parameter options (Placeholder)
Parameters can be freely changed within the main script.
There are several parameters that can be changed to include more stringent filtering for de novo peptides.


Path parameters 
|Parameter        |Default value| Description|
|-----------------|:-----------:|---------------|
|Default| True|                 If True, will look for setup file and overwrite manual filepaths.|  
|diamond_path| ..\CHEW\Setup\diamond\diamond.exe||
|diamond_folder| ..\CHEW\Setup\diamond\||
|ncbi_taxonomy_path|  ..\CHEW\Setup\ncbi_taxonomy\parsed_ncbi_taxonomy.tsv||        
|fasta_database_path  |||
|diamond_database_path|||

Path parameters specify which databases should be used. 


Filter parameters
|Parameter        |Default value| Description|
|-----------------|:-----------:|---------------|
|min_ALC_score|70                numeric, minimum required ALC score (Peaks score)|
|bit |25       |                   numeric, minimum required bitscore for alignment|
|freq_cut|5                      numeric, minimum lineage frequency filter for composition and DB creation
|Write_to_database| "Proteins"| do not make a database(False), use aligned proteins ("Proteins") use aligned taxids("Taxids").|
|DB_rank|"genus"|                 rank specificity  ("OX" "species" "genus" or "family") if "Taxids" is used for Write_to_database |
|Temporary_directory| ..\CHEW\|    folder for writing temporary DIAMOND indices |




#### Licensing

The pipeline is licensed with standard MIT-license. <br>
If you would like to use this pipeline in your research, please cite the following papers: 
      
- NovoLign: metaproteomic profiling by de novo sequence alignment <br>         

- Kleikamp, Hugo BC, et al. "Database-independent de novo metaproteomics of complex microbial communities." Cell Systems 12.5 (2021): 375-383.

- Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. "Fast and sensitive protein alignment using DIAMOND." Nature methods 12.1 (2015): 59-60.



#### Contact:
-Hugo Kleimamp (Developer): hugo.kleikamp@uantwerpen.be<br> 
-Martin Pabst (Co-Developer): M.Pabst@tudelft.nl<br>


#### Related repositories:
https://github.com/bbuchfink/diamond<br>
https://github.com/hbckleikamp/proteomic-database-prepper<br>
https://github.com/hbckleikamp/NCBI2Lineage<br>
https://github.com/hbckleikamp/GTDB2DIAMOND


