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


<br>
Using alignment parameters optimized to short peptide homology

LCAs optimized for short peptide homology

#### What does it do? (Placeholder)
<br>
Outputs




#### Running NovoLign (Placeholder)

<br>



#### What input files does it use? (Placenhlder)

<br>

#### What outputs does it generate? (Placeholder)

<br>

## Parameter options (Placeholder)
Parameters can be freely changed within the main script.
There are several parameters that can be changed to include more stringent filtering for de novo peptides.






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


