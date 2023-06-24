#!/usr/bin/env python

'''
make-chromosome-plasmid-table.py by Rohan Maddamsetti.

This script reads in ../results/best-prokaryotes.txt.
'''
import os

warning_count = 0
found_genomes = 0
accessions_not_found_in_refseq = []

with open("../results/chromosome-plasmid-table.csv",'w') as out_fh:
    header = "Organism,Strain,NCBI_Nucleotide_Accession,SequenceType,Annotation_Accession\n"
    out_fh.write(header)
    ## open the genome report file, and parse line by line.
    with open("../results/best-prokaryotes.txt", "r") as genome_report_fh:
        for i, line in enumerate(genome_report_fh):
            line = line.strip()
            if i == 0: ## get the names of the columns from the header.
                column_names_list = line.split('\t')
                continue ## don't process the header any further.
            fields = line.split('\t')
            organism = fields[0]
            strain = fields[-1]
            replicons = fields[8]
            ftp_path = fields[20]
            GBAnnotation = os.path.basename(ftp_path)
            my_annotation_file = "../results/gbk-annotation/" + GBAnnotation + "_genomic.gbff"
            ''' make sure that this file exists in the annotation directory--
            skip if this was not the case.
            this is important; we don't want to include genomes that were
            not in the search database in the first place. '''
            if not os.path.exists(my_annotation_file):
                warning_count += 1
                accessions_not_found_in_refseq.append(GBAnnotation)
                continue
            ## if we got here, then the RefSeq genbank annotation file was found.
            found_genomes += 1
            replicon_list = replicons.split(';')
            for s in replicon_list:
                s = s.strip()
                seq_id = s.split('/')[-1]
                if ':' in seq_id:
                    seq_id = seq_id.split(':')[-1]
                if s.startswith("chromosome"):
                    seq_type = "chromosome"
                elif s.startswith("plasmid"):
                    seq_type = "plasmid"
                else: ## only allow chromosomes and plasmids.
                    known_phage_seqs = ['CP003186.1','CP013974.1','CP019719.1',
                                        'GQ866233.1','CP018841.1','CP002495.1',
                                        'CP011970.1', 'CP027117.1', 'CP002889.1',
                                        'CP004084.1', 'CP011103.1', 'CP063968.1', 'CP062992.1']
                    if seq_id in known_phage_seqs:
                        continue ## pass silently.
                    else:
                        print("PROBABLE PHAGE (not chromosome or plasmid):")
                        print(s)
                        continue
                row_data = [organism, strain, seq_id, seq_type, GBAnnotation]
                ## replace all commas with semicolon to respect the csv output format.
                row_string = ','.join([x.replace(',',';') for x in row_data]) + '\n'
                out_fh.write(row_string)

## now print out warnings.
print("Warning: the following accessions were not found in RefSeq. ")
for GBAnnotation in accessions_not_found_in_refseq:
    print(GBAnnotation)
print("A total of " + str(warning_count) + " accessions were not found in RefSeq.")
print("The most likely explanation is that these missing accessions did not pass RefSeq Quality Control.")
print("A total of " + str(found_genomes) + " genomes were found in RefSeq and processed.")
