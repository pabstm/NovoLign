# NovoLign



This is the repository for the NovoLign pipeline, as described in:<br>
"NovoLign: metaproteomic profiling by de novo sequence alignment" 


The pipeline was established and tested with shotgun (meta)proteomics data obtained from Q Exactive Orbitrap Mass Spectrometers, using PEAKS generated de novo sequence lists. The generation of accurate de novo peptide sequence lists depends on high quality peptide sequencing spectra. NovoLign has been tested and developed in an Anaconda Spyder environment.

<img src="https://github.com/hbckleikamp/NovoLign/blob/main/images/workflow.svg" width="300" height="400" align="right">
<br>

#### What is NovoLign? 
NovoLign is a tool for rapid annotation of de novo sequenced peptides with homology alignment.
It uses DIAMOND for high-througput, error-tolerant annotation of de novo sequenced peptides.
As de novo sequencing is independent of database composition, it provides an unbiased alternative to conventional database searching methods. By aligning de novo peptides to general databases such as UniRef100, NovoLign can find related peptide sequences, which can aid in construction of sample specific databases, database and experiment quality control, and asses the completeness of references databases.


## Basic use

#### Setting up NovoLign 
NovoLign is designed to interface with multiple commonly used reference databases, belonging to NCBI, UniprotKB or GTDB.
Setup scripts are supplied in the folder `Setup` and can be run automatically with `run_Setup.py`
To setup protein databases for GTDB, utility scripts are supplied to download unzip and merge GTDB reference proteomes, and add taxonomies to fasta headers. For NCBI, a linearized taxonomy is constructed for specified taxonomic ranks using taxdump.
scripts are supplied for database curation, such as euqating I and L, and removing ambiguous amino acids.
Lastly, scripts are supplied for automated downloading of DIAMOND.

<br>

#### How does it work? 

The NovoLign pipeline consists of 5 parts:
1. *DIAMOND alignment* : Homology alignment is performed with parameters optimized for de novo sequencing errors. <br>
2. *Lowest Common Ancestor analysis (LCA)* : Different LCA algorithms are employed to maximize taxonomic specificity. <br>
3. *Taxonomy report* : Taxonomic quantification is performed with tabular and visual output. <br>
4. *Spectral quality report* : An assessment of experiment quality control is performed based on de novo scores and annotation rates.  <br>
5. *Database coverage report* : Optionally: to perform quality control on a reference database, outputs of de novo sequencing are compared against database searching outputs.  <br>




#### What does it do? 
<br>
After determining the best-performing parameter combinations for the NovoLign pipeline using synthetic communities, we evaluated its general practicability by studying a range of pure reference strains, enrichment cultures, and complex microbial communities. In addition to determining microbial composition using de novo sequence alignment, we analyzed for all experiments the fraction of high-quality fragmentation spectra that were matched during database searching. This provides an measure for the unmatched fraction during database searching, and therefore, for the completeness of the reference sequence databases used for database searching. De novo sequence alignment of the unmatched spectra moreover determines whether there are any taxonomies that are not covered by the reference sequence database used for database searching. Finally, the pipeline enables to construct a de novo focused UniRef100 reference sequence databases from the taxonomic composition determined by sequence alignment.![image](https://github.com/hbckleikamp/NovoLign/assets/49785660/9758d1e9-278c-4d99-975c-85cf3af66f8a)



#### Running NovoLign 
- NovoLign is designed as a single "tunable" python script.
- NovoLign does not offer command line options, but parameters can be altered in the main script.



#### What input files does it use? 
NovoLign is tested to work with .psm output formats from PEAKS de novo sequencing and DeepNovo.
Any tabular or .txt-like format can be supplied, provided it contains a column of peptide sequences with the header `Peptide`.
In default operation NovoLign will look for any folder starting with `Input_` within the NovoLign directory (see: Path parameters). Examples of input files are supplied in the folder `Input_p_yeast`
<br>

#### What outputs does it generate? 

<br>

## Parameter options (Placeholder)
Parameters can be freely changed within the main script.
There are several parameters that can be changed to include more stringent filtering for de novo peptides.


Path parameters specify which databases should be used. 
|Parameter        |Default value| Description|
|-----------------|:-----------:|---------------|
|Default| True|                 If True, will look for setup file and overwrite manual filepaths.|  
|input_files| ..\CHEW\*Input_*| Location of input folder
|diamond_path| ..\CHEW\Setup\diamond\diamond.exe |Location of DIAMOND executable|
|diamond_folder| ..\CHEW\Setup\diamond\ |Location of DIAMOND folder|
|Temporary_directory| ..\CHEW\ |    Folder for writing temporary DIAMOND indices |
|ncbi_taxonomy_path|  ..\CHEW\Setup\ncbi_taxonomy\parsed_ncbi_taxonomy.tsv|Location of linear NCBI taxonomy|        
|fasta_database_path  ||Location of database fasta file|
|diamond_database_path||Location of database dmnd file|


Other parameters specify cutoffs for de novo score, alignment score and taxon frequency, as how use NovoLign to perform database construction.
|Parameter        |Default value| Description|
|-----------------|:-----------:|---------------|
|min_ALC_score| 70 |               numeric, minimum required ALC score (Peaks score)|
|bit | 25 |                   numeric, minimum required bitscore for alignment|
|freq_cut|5 |                     numeric, minimum lineage frequency filter for composition and DB creation|
|Write_to_database| "Proteins"| do not make a database(False), use aligned proteins ("Proteins") use aligned taxids("Taxids").|
|DB_rank|"genus"|                 rank specificity  ("OX" "species" "genus" or "family") if "Taxids" is used for Write_to_database |





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


