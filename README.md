# AR-gene-plasmid-analysis by Rohan Maddamsetti and Vincent Huang

## Python requirements: Python 3.6+, biopython, tqdm

## External software requirements:
Mash 2.3: https://github.com/marbl/Mash  
kallisto 0.46: https://pachterlab.github.io/kallisto/about  
DIAMOND 2.1.6: http://www.diamondsearch.org  
Assembly Dereplicator 0.3.1: https://github.com/rrwick/Assembly-Dereplicator  

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

sbatch --mem=16G -t 24:00:00 --wrap="python fetch-gbk-annotation.py"  
sbatch --mem=16G -t 24:00:00 --wrap="python fetch-assembly-stats.py"  

Now run the following scripts on DCC. Some run
quite quickly, so no need to submit them to a partition on DCC--
just run them in an interactive session on DCC.

python make-chromosome-plasmid-table.py  
python make-gbk-annotation-table.py ## this runs for ~35 min on DCC.

## double-check assembly quality on DCC.  
python run-QC-and-make-assembly-stats-table.py  

## this runs for ~8h on DCC.
sbatch --mem=16G -t 24:00:00 --wrap="python count-cds.py"  

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
prokaryotes-with-plasmids.txt  
genome-assembly-metadata.csv  


Locally, download fasta sequences for all genomes, and make a list of dereplicated
sequences. This runs overnight, and uses a lot of memory (100Gb!):
python fetch-and-dereplicate-seqs.py

Then, run the follow scripts to annotate the genomes, and to cross-check
the computational annotation against a subset of annotations that were conducted manually.  

python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv  

python check-ecological-annotation.py  

### CARD and mobileOG-db analyses.

After downloading duplicate-proteins.csv and all-proteins.csv, run the following locally:  
python protein_csv_to_fasta.py  
python protein_csv_to_fasta.py --ignore-singletons  

Then download the CARD database locally into a folder called "../data/card-data", relative to the directory
containing this source code, and download the mobileOG-db database into a folder called
"../data/mobileOG-db_beatrix-1-6_v1_all" relative to the directory containing this source code.

Make sure all the paths in this next script make sense, and run the following:  
python search-CARD-and-mobileOG-db.py  

Then run:  
python parse-DIAMOND-results.py.  

### Now, all the data analysis in ARG-duplications.R should run.