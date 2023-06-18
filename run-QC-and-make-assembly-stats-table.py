#!/usr/bin/env python

'''
run-QC-and-make-assembly-stats-table.py by Rohan Maddamsetti

Print out a table of assembly metadata, subsettting on assemblies
that are deposited in RefSeq (i.e. passed RefSeq QC and no assembly anomalies),
and are complete genomes (no spanned or unspanned gaps, and no unplaced scaffolds)

'''

import os
from tqdm import tqdm

assembly_stats_dir = "../results/refseq-assembly-statistics/"

with open("../results/genome-assembly-metadata.csv","w") as out_fh:
    header = "Annotation_Accession,Organism_name,Assembly_level,Assembly_method,Genome_coverage,Sequencing_technology,Genbank_assembly_accession,Refseq_assembly_accession,total_length,spanned_gaps,unspanned_gaps,total_gap_length,scaffold_count\n"
    out_fh.write(header)

    stats_files = [x for x in os.listdir(assembly_stats_dir) if x.endswith("_assembly_stats.txt")]
    
    for x in tqdm(stats_files):
        assembly_stats_path = assembly_stats_dir + x
        annotation_accession = x.split("_assembly_stats.txt")[0]
        organism_name = "NA"
        assembly_level = "NA"
        assembly_method = "NA"
        genome_coverage = "NA"
        sequencing_technology = "NA"
        genbank_accession = "NA"
        refseq_accession = "NA"
        ## key statistics on the assembly for QC.
        total_length = "NA"
        spanned_gaps = "NA"
        unspanned_gaps = "NA"
        total_gap_length = "NA"
        scaffold_count = "NA"

        QC_failed = False
        
        with open(assembly_stats_path,'rt') as gbk_fh:
            for line in gbk_fh:
                line = line.strip()
                if "Excluded from RefSeq" in line:
                    QC_failed = True
                    break
                elif "Assembly anomaly" in line:
                    QC_failed = True
                    break
                elif "unlocalized-scaffold" in line:
                    if "total-length" in line:
                        unlocalized_scaffold_length = int(line.split()[-1])
                        if unlocalized_scaffold_length != 0:
                            QC_failed = True
                            break
                elif line.startswith("# Organism name:"):
                    organism_name = line.split("# Organism name:")[-1].strip()
                elif line.startswith("# Assembly level:"):
                    if "Complete Genome" not in line:
                        QC_failed = True
                        break
                    else:
                        assembly_level = "Complete Genome"
                elif line.startswith("# Assembly method:"):
                    assembly_method = line.split("# Assembly method:")[-1].strip()
                elif line.startswith("# Genome coverage:"):
                    genome_coverage = line.split("# Genome coverage:")[-1].strip()
                elif line.startswith("# Sequencing technology:"):
                    sequencing_technology = line.split("# Sequencing technology:")[-1].strip()
                elif line.startswith("# GenBank assembly accession:"):
                    genbank_accession = line.split("# GenBank assembly accession:")[-1].strip()
                elif line.startswith("# RefSeq assembly accession:"):
                    refseq_accession = line.split("# RefSeq assembly accession:")[-1].strip()
                elif line.startswith("all"):
                    if "total-length" in line:
                        total_length = line.split()[-1]
                    ## need the leading tab '\t' since 'spanned-gaps' is a substring of 'unspanned-gaps'
                    elif "\tspanned-gaps" in line:
                        spanned_gaps = line.split()[-1]
                        if int(spanned_gaps) != 0:
                            QC_failed = True
                            break
                    ## ditto
                    elif "\tunspanned-gaps" in line:
                        unspanned_gaps = line.split()[-1]
                        if int(unspanned_gaps) != 0:
                            QC_failed = True
                            break
                    elif "total-gap-length" in line:
                        total_gap_length = line.split()[-1]
                        if int(total_gap_length) != 0:
                            QC_failed = True
                            break
                    elif "scaffold-count" in line:
                        scaffold_count = line.split()[-1]
                else:
                    continue
            if not QC_failed:
                ## remove all commas found in each entry
                rowdata = [x.replace(",", "") for x in [annotation_accession, organism_name, assembly_level, assembly_method, genome_coverage, sequencing_technology, genbank_accession, refseq_accession, total_length, spanned_gaps, unspanned_gaps, total_gap_length, scaffold_count]]
                row = ','.join(rowdata) + '\n'
                out_fh.write(row)
