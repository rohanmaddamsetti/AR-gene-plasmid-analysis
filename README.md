# AR-gene-plasmid-analysis by Rohan Maddamsetti and Vincent Huang

## Python requirements: Python 3.6, biopython, tqdm 


First, download plasmids.txt and prokaryotes.txt into
../data/GENOME_REPORTS:

wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

Then, filter the prokaryotes genome data for those that have either complete genomes or both
plasmids and chromosomes:

python filter-genome-reports.py > ../results/best-prokaryotes.txt

Then, fetch genome annotation for each row in best-prokaryotes.txt,
and fetch the protein-coding genes for all chromosomes and plasmids for
each row in best-prokaryotes.txt.

Both steps can be done at the same time on the Duke Compute Cluster (DCC).
And make sure these scripts are called from the src directory.
fetch-gbk-annotation and fetch-genome-and-plasmid-cds.py run overnight..

sbatch --mem=16G -t 24:00:00 --wrap="python fetch-gbk-annotation.py"  
sbatch --mem=16G -t 24:00:00 --wrap="python fetch-genome-and-plasmid-cds.py"  

Once the data has downloaded, run the following scripts. Some run
quite quickly, so no need to submit them to a partition on DCC--
just run them in an interactive session on DCC.

python make-chromosome-plasmid-table.py  
python make-gbk-annotation-table.py ## this runs for ~35 min on DCC.  

## this runs for 20 min on DCC. 
python count-cds.py > ../results/protein_db_CDS_counts.csv

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
protein_db.faa  
prokaryotes-with-plasmids.txt  

Then, run the follow scripts to annotate the genomes, and to cross-check
the computational annotation against a subset of annotations that were conducted manually.

python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv  

python check-ecological-annotation.py  

