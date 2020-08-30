#!/usr/bin/env python

'''
make-chromosome-plasmid-table.py by Rohan Maddamsetti.

This script reads in ../results/prokaryotes-with-plasmids.txt.
'''
import os

with open("../results/chromosome-plasmid-table.csv",'w') as out_fh:
    header = "Organism,Strain,NCBI_Nucleotide_Accession,SequenceType,Annotation_Accession\n"
    out_fh.write(header)
    ## open the genome report file, and parse line by line.
    with open("../results/prokaryotes-with-plasmids.txt", "r") as genome_report_fh:
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
            my_annotation_file = "../results/gbk_annotation/" + GBAnnotation + "_truncated.gbff"
            ''' make sure that this file exists in the annotation directory--
            skip if this was not the case.
            this is important; we don't want to include genomes that were
            not in the search database in the first place. '''
            if not os.path.exists(my_annotation_file):
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
                    print('WEIRD CASE:')
                    print(s)
                    continue
                row_string = ','.join([organism, strain, seq_id, seq_type, GBAnnotation]) + '\n'
                out_fh.write(row_string)
