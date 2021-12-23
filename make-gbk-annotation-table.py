#!/usr/bin/env python

'''
make-gbk-annotation-table.py by Rohan Maddamsetti
'''

import os
import gzip
from tqdm import tqdm

gbk_annotation_dir = "../results/gbk-annotation/"

with open("../results/gbk-annotation-table.csv","w") as out_fh:
    header = "Annotation_Accession,host,isolation_source\n"
    out_fh.write(header)

    gbk_files = [x for x in os.listdir(gbk_annotation_dir) if x.endswith("_genomic.gbff.gz")]
    
    for x in tqdm(gbk_files):
        gbk_path = gbk_annotation_dir + x
        annotation_accession = x.split("_genomic.gbff.gz")[0]
        with gzip.open(gbk_path,'rt') as gbk_fh:
            host = "NA"
            isolation_source = "NA"
            ## use a buffer to store the line, to handle cases where
            ## the annotation spans multiple lines.
            in_host_field = False
            in_isolation_source_field = False
            line_buffer = []
            for line in gbk_fh:
                line = line.strip()
                ## We're going to delete all double-quote characters,
                ## and replace all commas with semicolons so that they
                ## don't clash with the csv format.
                ## BUT: make sure that the line terminates with a double-quote--
                ## otherwise, we need to slurp up data from multiple lines.
                if line.startswith("/host"):
                    line_annot = line.split('=')[-1].replace('\"','').replace(',',';')
                    if line.endswith('"'):
                        host = line_annot
                    else: ## have to look at the next line too.
                        line_buffer.append(line_annot)
                        in_host_field = True
                elif in_host_field and len(line_buffer):
                    line_annot = line.replace('\"','').replace(',',';')
                    line_buffer.append(line_annot)
                    if line.endswith('"'): ## flush the line buffer and reset the flag.
                        host = ' '.join(line_buffer)
                        line_buffer = []
                        in_host_field = False
                elif line.startswith("/isolation_source"):
                    line_annot = line.split('=')[-1].replace('\"','').replace(',',';')
                    if line.endswith('"'):
                        isolation_source = line_annot
                    else: ## then have to look at next line too.
                        line_buffer.append(line_annot)
                        in_isolation_source_field = True
                elif in_isolation_source_field and len(line_buffer):
                    line_annot = line.replace('\"','').replace(',',';')
                    line_buffer.append(line_annot)
                    if line.endswith('"'): ## flush the line buffer and reset flag.
                        isolation_source = ' '.join(line_buffer)
                        line_buffer = []
                        in_isolation_source_field = False
                ## break if we have the annotation.
                if (host != "NA") and (isolation_source != "NA"): break
                ## also break if we're looking at gene annotation since
                ## there's insufficient isolate annotation in this file.
                if line.startswith("/gene"): break ## NOTE: has to be '/gene'
                ## because we don't want to match the Genome Annotation Data in the
                ## COMMENT field of the Genbank metadata.
                if line.startswith("ORIGIN"): break ## break if we got to the first fasta seq
            row = ','.join([annotation_accession, host, isolation_source]) + '\n'
            out_fh.write(row)
