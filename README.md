# NovoLign



This is the repository for the NovoLign pipeline, as described in:<br>
"Kleikamp, Hugo BC, et al. "NovoLign: metaproteomics by sequence alignment." ISME communications 4.1 (2024): ycae121." 

This github page will be continuously updated. To use the original NovoLign version as was published, 
please refer to the archived github page:
https://github.com/hbckleikamp/NovoLign/


The pipeline was established and tested with shotgun proteomics and metaproteomics data obtained from different projectss. The measurments were performed using Orbitrap mass spectrometers, de novo sequence lists were generated using PEAKS Studio. The generation of accurate de novo peptide sequence lists depends on high quality peptide sequencing spectra. NovoLign has been tested and developed using the Anaconda Spyder environment.
It is highly recommended to use Python package versions as specified in the "NovoLign_requirements_17032024" text file to ensure smooth operation of the pipeline. Newer versions may lead to inconsistent outcomes and unexpected errors.

#### What is NovoLign? 
NovoLign is a tool for rapid annotation of de novo sequenced peptides from complete metaproteomics experiments uisng homology alignment. It uses DIAMOND for high-througput, error-tolerant annotation of de novo sequenced peptides.
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
4. *Spectral quality report* : Optionally: An assessment of experiment quality control is performed based on de novo scores and annotation rates.  <br>
5. *Database coverage report* : Optionally: to perform quality control on a reference database, outputs of de novo sequencing are compared against database searching outputs.  <br>




#### What does it do? 
After determining the best-performing parameter combinations for the NovoLign pipeline using synthetic communities, we evaluated its general practicability by studying a range of pure reference strains, enrichment cultures, and complex microbial communities. In addition to determining microbial composition using de novo sequence alignment, we analyzed for all experiments the fraction of high-quality fragmentation spectra that were matched during database searching. This provides an measure for the unmatched fraction during database searching, and therefore, for the completeness of the reference sequence databases used for database searching. De novo sequence alignment of the unmatched spectra moreover determines whether there are any taxonomies that are not covered by the reference sequence database used for database searching. Finally, the pipeline enables to construct a de novo focused UniRef100 reference sequence databases from the taxonomic composition determined by sequence alignment.


#### Running NovoLign ####
NovoLign can be run in two ways:<br><br>
1.Executing Novolign from command line: 
Here arguments are supplied in the command line.
<br><br>
2.Executing NovoLign directly from the most recent script (NovoLign_workflow_v27112024.py):
Here filepaths are changed manually each time in the script.

NovoLign is tested to work with .psm output formats from PEAKS de novo sequencing and DeepNovo.
Any tabular or .txt-like format can be supplied, provided it contains a column of peptide sequences with the header `Peptide`.
If an "input folder" is supplied , its contents should adhere to the following structure (See Example: Input_p_Yeast)
- de_novo_file=..\\*de novo* *peptides.csv
- database_searching_file=..\\*psm.csv
- fasta_database=..\\*.fasta (used in database qc only)
- taxa=..\\*TaxIDs.txt (used in database qc only)

### Execute Novolign from command line ###

Example commands, executed from the NovoLign directory after running Setup.py:

Base NovoLign usage (just alignment and LCA)
```
python "NovoLign_CLI.py" -i "Input_p_Yeast/de novo peptides.csv" -d "Setup/Swiss-Prot/uniprot_sprot_NoAmb_IJeqL.dmnd"
```

With database construction:
```
python "NovoLign_CLI.py" -i "Input_p_Yeast/de novo peptides.csv" -d "Setup/Swiss-Prot/uniprot_sprot_NoAmb_IJeqL.dmnd" -f "Setup/Swiss-Prot/uniprot_sprot_NoAmb_IJeqL.fasta"
```

With comparison to database searching output:
```
python "NovoLign_CLI.py" -i "Input_p_Yeast/de novo peptides.csv" -di "Input_p_Yeast/DB search psm.csv" -d "Setup/Swiss-Prot/uniprot_sprot_NoAmb_IJeqL.dmnd" 
```


<br>**Main arguments**
|Argument      |Default value| Description|
|-----------------|:-----------:|---------------|
|-i| n.a. | Required, --input_file: one or more comma separated filepaths (or an input folder with NovoLign folder structure) |
|-d| n.a. | Required, --diamond_database_path: Path to DIAMOND database used for alignment |
|-o| n.a. | Optional, --Output_directory: Path to output folder|

<br>**Performance arguments**
|Argument      |Default value| Description|
|-----------------|:-----------:|---------------|
|-ALC| 70 | Optional, --min_ALC_score: Minimum ALC score (PEAKS specific score)|
|-bit| 25 | Optional, -min_bit_score: Minimum bitscore (DIAMOND alignment score)|
|-freq| 5 | Optional, --freq_cut: Minium taxa frequency for denoising|
|-lcas| ['lca','bitlca','wcla'] | Optional, which lca algorithms to use: lca (conventional lca), bitlca (bistcore weighted lca), wlca (weighted lca)|

<br>**Database construction**
|Argument      |Default value| Description|
|-----------------|:-----------:|---------------|
|-f| n.a. | Optional, --fasta_database_path: Path to fasta database (required for database construction) |
|-DBwrite| "Proteins" | Optional, Used only when -f is defined for database construction. Options: (False, 'Proteins','Taxids'): do not make a database(False), use aligned proteins ('Proteins') use aligned taxids('Taxids').|
|-DB_rank| "genus" | Optional, Used only in database construction with -DBwrite 'Taxids', selects the taxonomic rank for database construction, Options: 'OX' 'species' 'genus' or 'family' |

<br>**Database searching comparison**
|Argument      |Default value| Description|
|-----------------|:-----------:|---------------|
|-di| n.a. | Optional, --database_searching_file: database searching file, for comparison with de novo sequenced peptides, should be tabular and contain the column 'Peptide' |
|-taxa| n.a. | Optional, Specified taxa for the visual comparison of database searching and de novo sequencing outputs |
|-o| n.a. | Optional, --Output_directory: Path to output folder|

<br>**Default filepaths**
|Argument      |Default value| Description|
|-----------------|:-----------:|---------------|
|diamond_path| ..\NovoLign\Setup\diamond\diamond.exe |Location of DIAMOND executable|
|diamond_folder| ..\NovoLign\Setup\diamond\ |Location of DIAMOND folder|
|Temporary_directory| ..\NovoLign\ |    Folder for writing temporary DIAMOND indices |
|ncbi_taxonomy_path|  ..\NovoLign\Setup\ncbi_taxonomy\parsed_ncbi_taxonomy.tsv|Location of linear NCBI taxonomy|        







### Execute Novolign by running script ###
This requires manually editing filepaths to your input folder and diamond database in recent script (NovoLign_workflow_v27112024.py).
Examples of input files are supplied in the folder `Input_p_yeast`
<br>

#### What outputs does it generate? 

NovoLign generates several output files divided over different folders.
|Folder           | Section     |       Contents|
|-----------------|:-----------:|---------------|
|diamond_fasta| 1 | generated fasta file from input peptide sequences file for DIAMOND alignment|
|diamond_alignments| 1 | DIAMOND alignment|
|lca| 2 | different LCA outputs: conventional (CON) weighted (W), bitscore weighted (BIT) |
|composition| 3 | taxonomic composition of input sample (for different LCAs) including level of decoy matches|
|experiment_qc| 4 | comparison of spectral annotation rates of input sample by NovoLign and database searching, for spectra at different quality levels|
|database_qc| 5 | comparison of taxonomic composition of input sample obtained by NovoLign to taxonomic composition obtained by database searching |
|psms| 5 | Final PSMs format output with NovoLign annotation |

<br>
Example output vizualisation for experiment quality control, which compares annotation rates for de novo sequenced PSMs at different de novo scores (ALC%) compared to database searching, and checks if the aligned peptide sequences are the same as the peptide sequences found in database searching.
<br>
<br>
<p align="center">
      <img src="https://github.com/hbckleikamp/NovoLign/blob/main/images/DB search psm_bins.png" width="45%" height="300" align="left">
&nbsp; &nbsp; &nbsp; &nbsp;
      <img src="https://github.com/hbckleikamp/NovoLign/blob/main/images/de novo peptides_bins.png" width="45%" height="300" align="right">
</p>
<br clear="left"/>

<br>
Example output vizualisation for database quality control, which compares the taxonomic composition of PSMs found exclusively in de novo sequencing (DN_only), all de novo PSMs (DN_all), all database searching PSMs (DB_all) and all PSMs unique to low scoring database searching hits (DBLQ_only).
<br>
<br>
<p align="center">
      <img src="https://github.com/hbckleikamp/NovoLign/blob/main/images/DB_vs_DN_genusabsolute_topX.png" width="45%" height="300" align="left">
&nbsp; &nbsp; &nbsp; &nbsp;
      <img src="https://github.com/hbckleikamp/NovoLign/blob/main/images/DB_vs_DN_genusnormalized_topX.png" width="45%" height="300" align="right">
</p>
<br clear="left"/>




#### Licensing

The pipeline is licensed with standard MIT-license. <br>
If you would like to use this pipeline in your research, please cite the following papers: 
      
- Kleikamp, Hugo BC, et al. "NovoLign: metaproteomics by sequence alignment." ISME communications 4.1 (2024): ycae121..<br>         

- Kleikamp, Hugo BC, et al. "Database-independent de novo metaproteomics of complex microbial communities." Cell Systems 12.5 (2021): 375-383.

- Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. "Fast and sensitive protein alignment using DIAMOND." Nature methods 12.1 (2015): 59-60.



#### Contact:
-Hugo Kleimamp (Developer): hugo.kleikamp@uantwerpen.be<br> 
-Martin Pabst (Co-Developer): m.pabst@tudelft.nl<br>


#### Related repositories:
https://github.com/bbuchfink/diamond<br>
https://github.com/hbckleikamp/proteomic-database-prepper<br>
https://github.com/hbckleikamp/NCBI2Lineage<br>
https://github.com/hbckleikamp/GTDB2DIAMOND<br>
https://github.com/hbckleikamp/De-Novo-ErrorSIM


