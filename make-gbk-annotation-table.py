#!/usr/bin/env python

'''
make-gbk-annotation-table.py by Rohan Maddamsetti
'''

import os
import gzip

gbk_annotation_dir = "../results/gbk-annotation/"

with open("../results/gbk-annotation-table.csv","w") as out_fh:
    header = "Annotation_Accession,host,isolation_source\n"
    out_fh.write(header)

    gbk_files = [x for x in os.listdir(gbk_annotation_dir) if x.endswith("_genomic.gbff.gz")]
    for x in gbk_files:
        gbk_path = gbk_annotation_dir + x
        annotation_accession = x.split("_genomic.gbff.gz")[0]
        with gzip.open(gbk_path,'rt') as gbk_fh:
            host = "NA"
            isolation_source = "NA"
            for line in gbk_fh:
                line = line.strip()
                ## We're going to delete all double-quote characters,
                ## and replace all commas with semicolons so that they
                ## don't clash with the csv format.
                if line.startswith("/host"):
                    host = line.split('=')[-1].replace('\"','').replace(',',';')
                elif line.startswith("/isolation_source"):
                    isolation_source = line.split('=')[-1].replace('\"','').replace(',',';')
                if (host != "NA") and (isolation_source != "NA"): break
            row = ','.join([annotation_accession, host, isolation_source]) + '\n'
            out_fh.write(row)
