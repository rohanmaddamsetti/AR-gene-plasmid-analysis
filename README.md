# AR-gene-plasmid-analysis by Rohan Maddamsetti and Vincent Huang

## Software requirements
## Python 3.6+, biopython, tqdm  
## R 4.0  
## Julia 1.8.5  
### Mash 2.3: https://github.com/marbl/Mash  
### Assembly Dereplicator 0.3.1: https://github.com/rrwick/Assembly-Dereplicator  
### DIAMOND 2.1.6: http://www.diamondsearch.org  
### kallisto 0.48.0: https://pachterlab.github.io/kallisto/  

Make a top-level directory with three directories inside, named "data", "results", and "src".  
Now copy all source code files in this repository into "src".  

Now, download prokaryotes.txt into ../data/GENOME_REPORTS:  

wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt  

Then, filter the prokaryotes.txt genome data for those that have complete genomes,
and replace "GCA" with "GCF" throughout this file, so that RefSeq data and not Genbank data
is accessed in all downstream steps:  

python filter-genome-reports.py > ../results/best-prokaryotes.txt  

Then, fetch genome annotation for each row in best-prokaryotes.txt,
fetch the protein-coding genes for all chromosomes and plasmids for
each row in best-prokaryotes.txt,
and fetch the assembly statistics for quality control.
We want to analyze genomes with high-quality, complete genome assemblies.  

These steps can be done at the same time on the Duke Compute Cluster (DCC).
And make sure these scripts are called from the src directory.
fetch-gbk-annotation runs for several hours.  

First, make some output directories.  
mkdir ../results/gbk-annotation/  
mkdir ../results/refseq-assembly-statistics/  

Then:  

sbatch --mem=16G -t 36:00:00 --wrap="python fetch-gbk-annotation.py"  
sbatch --mem=16G -t 36:00:00 --wrap="python fetch-assembly-stats.py"  

Now run the following scripts on DCC. Some run
quite quickly, so no need to submit them to a partition on DCC--
just run them in an interactive session on DCC.

python make-chromosome-plasmid-table.py  
python make-gbk-annotation-table.py ## this runs for ~35 min on DCC.

## double-check assembly quality on DCC.  
python run-QC-and-make-assembly-stats-table.py  

## this runs for ~12h on DCC.
sbatch --mem=16G -t 24:00:00 --wrap="python count-cds.py"  
sbatch --mem=16G -t 24:00:00 --wrap="python count-proteins-and-replicon-lengths.py"  

## this runs for ~36h on DCC.
sbatch --mem=16G -t 48:00:00 --wrap="python tabulate-proteins.py"  

## this runs for ~36h on DCC.
sbatch --mem=16G -t 48:00:00 --wrap="python tabulate-proteins.py --ignore-singletons"  

Then, copy the following files from the results/
directory onto my local machine (same directory name and file structure).

duplicate-proteins.csv  
all-proteins.csv  
protein_db_CDS_counts.csv  
gbk-annotation-table.csv  
chromosome-plasmid-table.csv  
genome-assembly-metadata.csv  
replicon-lengths-and-protein-counts.csv


Locally, download fasta sequences for all genomes, and make a list of dereplicated
sequences. This runs overnight, and uses a lot of memory (100Gb!):  
python fetch-and-dereplicate-seqs.py

Then, run the follow scripts to annotate the genomes, and to cross-check
the computational annotation against a subset of annotations that were conducted manually.  

python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv  

python check-ecological-annotation.py  

Finally, run the following script to find any genomes with chromosomes smaller than
a plasmid in that genome:  

python find-bad-replicons.py > ../results/bad-replicons.csv  


### CARD and mobileOG-db analyses.

After downloading duplicate-proteins.csv and all-proteins.csv, run the following locally:  
python protein_csv_to_fasta.py  
python protein_csv_to_fasta.py --ignore-singletons  

Then download the CARD database locally into a folder called "../data/card-data", relative to the directory
containing this source code, and download the mobileOG-db database into a folder called
"../data/mobileOG-db_beatrix-1-6_v1_all" relative to the directory containing this source code.

These databases can be downloaded from the following links:
https://card.mcmaster.ca/download  
https://mobileogdb.flsi.cloud.vt.edu/entries/database_download  

Make sure all the paths in this next script make sense, and run the following:  
python search-CARD-and-mobileOG-db.py  

Then run:  
python parse-DIAMOND-results.py.  


### linkage analysis of ARGs and MGEs.

Run the following script to count how often ARGs and MGEs are physical
neighbors (i.e. show genetic linkage) in these microbial genomes:

python count-ARG-MGE-adjacencies.py

I ran this on DCC as follows:
sbatch -p scavenger -t 24:00:00 --mem=4G --wrap="python count-ARG-MGE-adjacencies.py"  

The count-ARG-MGE-adjacencies.py script makes the following file:  
../results/ARG-MGE-adjacency-counts.csv  

Now run join-duplications.py to find larger regions of duplicated genes within each genome in the dataset.  
python join-duplications.py

This makes the following file:  
../results/joined-duplicate-proteins.csv  


### analysis of clinical antibiotic-resistance isolates.

Run the following script:
python tabulate-proteins-in-clinical-genomes.py  

This script makes the following files:  
../results/Duke-ESBL-all-proteins.csv  
../results/Duke-ESBL-duplicate-proteins.csv  
../results/BARNARDS-all-proteins.csv  
../results/BARNARDS-duplicate-proteins.csv  
../results/Mahmud2022-all-proteins.csv  
../results/Mahmud2022-duplicate-proteins.csv  
../results/Hawkey2022-all-proteins.csv  
../results/Hawkey2022-duplicate-proteins.csv  

Now run the following script to analyze plasmid copy number in the
dataset from Hawkey et al. (2022).  

python plasmid-copy-number-pipeline.py  

This script requires kallisto to be available from the command-line.  
It generates the following files used by ARG-duplications.R:  
../results/Hawkey2022_ARG_copy_numbers.csv  
../results/Hawkey2022_chromosome_plasmid_copy_numbers.csv  

### Finally -- now, all the data analysis in ARG-duplications.R should run.

### Analysis of transposases associated with duplicated ARGs.  

The julia script cluster-transposases.jl requires two input files that are generated by
ARG-duplications.R.
These files are:  
../results/Ecoli-transposases-in-dup-regions-with-ARGs.csv  
../results/transposases-in-dup-regions-with-ARGs.csv  

ARG-duplications.R also requires an output file generated by cluster-transposases.jl:  
../results/merged_transposases-in-dup-regions-with-ARGs.csv  

To avoid circular dependencies, ARG-duplications.R makes a system call to
run "julia cluster-transposases.jl".

Once ARG-duplications.R has run successfully, cluster-transposases.jl can also
be run as a stand-along script (since now its inputs exist on disk).

