#!/usr/bin/env python

'''
make-chromosome-plasmid-table.py by Rohan Maddamsetti.

This script reads in ../results/best-prokaryotes.txt.
'''
import os

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
            my_annotation_file = "../results/gbk-annotation/" + GBAnnotation + "_genomic.gbff.gz"
            ''' make sure that this file exists in the annotation directory--
            skip if this was not the case.
            this is important; we don't want to include genomes that were
            not in the search database in the first place. '''
            if not os.path.exists(my_annotation_file):
                print("ERROR:")
                print(GBAnnotation + " NOT FOUND:")
                print("searched for: " + my_annotation_file)
                print()
                continue
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
                                        'CP011970.1']
                    if seq_id in known_phage_seqs:
                        continue ## pass silently.
                    else:
                        print("WEIRD CASE:")
                        print(s)
                        continue
                row_data = [organism, strain, seq_id, seq_type, GBAnnotation]
                ## replace all commas with semicolon to respect the csv output format.
                row_string = ','.join([x.replace(',',';') for x in row_data]) + '\n'
                out_fh.write(row_string)
