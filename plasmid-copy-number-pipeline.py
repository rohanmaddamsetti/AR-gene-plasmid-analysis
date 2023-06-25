#!/usr/bin/env python

"""
plasmid-copy-number-pipeline.py by Rohan Maddamsetti.

For each genome in the Hawkey et al. 2022 dataset:

1) Download sequencing Illumina reads.
2) generate a fasta file of gene sequences from the genbank annotation.
     mark each gene by chromosome or plasmid, and give the replicon an ID.
     Also annotate the product field, so that ARGs and MGE-associated genes
     can be scored.
3) Run Kallisto Index on the fasta file of gene sequences.
4) Run kallisto quant on the index and sequencing reads.
5) estimate copy number of each replicon by averaging over genes,
    and make a table of copy number for all plasmids and chromosomes.
6) estimate copy number of all ARGs in the genome, relative to chromosome.
     Specifically, divide est_counts by length for each ARG to get coverage per bp,
     and divide by chromosome coverage per bp to get 
     copy number relative to chromosome.
    Then make a table of ARG copy number relative to chromosome.
7) Make a table of chromsome and plasmid lengths, based on
    the reference genomes.

NOTE: GCF_026154285.1_ASM2615428v1 did not have any reads pseudoalign.
I checked, and I did download the right SRA reads file. But it states at the following URL
that this RefSeq ID is suppressed. So maybe something wrong with this one?
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026154285.1/
"""

import subprocess
import os
import gzip
import re
from Bio import SeqIO


def download_fastq_reads(SRA_data_dir, SRA_accession_list_file):
    """
    use sra-toolkit to download the SRA accession data in
    SRA_accession_list_file into SRA_data_dir.
    """
    sra_numbers = []
    with open(SRA_accession_list_file, "r") as SRA_acc_fh:
        for SRA_acc in SRA_acc_fh:
            sra_numbers.append(SRA_acc.strip())
            
    for sra_id in sra_numbers:
        """
        the sra_id has to be the last part of the directory.
        see documentation here:
        https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
        """
        sra_dir_path = os.path.join(SRA_data_dir, sra_id)
        if os.path.exists(sra_dir_path): continue
        prefetch_args = ["prefetch", sra_id, "-O", sra_dir_path]
        print (" ".join(prefetch_args))
        subprocess.run(prefetch_args)

    ## have to change working directory for fasterq_dump
    my_cwd = os.getcwd()
    os.chdir(SRA_data_dir)
    for sra_id in sra_numbers:
        ## paired-end files for kallisto.
        sra_fastq_file_1 = sra_id + "_1.fastq"
        sra_fastq_file_2 = sra_id + "_2.fastq"
        if os.path.exists(sra_fastq_file_1) and os.path.exists(sra_fastq_file_2):
            continue
        else:
            print ("Generating fastq for: " + sra_id)
            ## run with 10+ threads  (default is 6 threads).
            fasterq_dump_args = ["fasterq-dump", "--threads", "10", sra_id]
            print(" ".join(fasterq_dump_args))
            subprocess.run(fasterq_dump_args)
    ## now change back to original working directory.
    os.chdir(my_cwd)
    return


def generate_fasta_reference_for_kallisto(gbk_gz_path, outfile):
    """ NOTE: this code depends on some assumptions that are true
    of the Hawkey et al. 2022 (by manual verification) that may not
    generalize for datasets uploaded by other people.
    """
    with open(outfile, "w") as outfh:
        with gzip.open(gbk_gz_path,'rt') as gbk_fh:
            SeqID = None
            SeqType = None
            for record in SeqIO.parse(gbk_fh, "genbank"):
                SeqID = record.id
                if "complete" in record.description:
                    if "plasmid" in record.description:
                        SeqType = "plasmid"
                    elif "chromosome" in record.description:
                        SeqType = "chromosome"
                    else:
                        continue
                else:
                    continue
                for feature in record.features:
                    ## only analyze protein-coding genes.
                    if feature.type != "CDS": continue
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    ## Important: for kallisto, we need to replace spaces with underscores in the product annotation field.
                    product = feature.qualifiers["product"][0].replace(" ","_")
                    DNAseq = feature.extract(record.seq)
                    header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"locus_tag="+locus_tag,"product="+product])
                    outfh.write(header + "\n")
                    outfh.write(str(DNAseq) + "\n")
    return


def make_Hawkey_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_outdir):    
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        print("making: ", fasta_outfile)
        generate_fasta_reference_for_kallisto(gzpath, fasta_outfile)
    return


def make_Hawkey_kallisto_indices(kallisto_ref_dir, kallisto_index_dir):
    ref_fasta_filelist = [x for x in os.listdir(kallisto_ref_dir) if x.endswith(".fna")]
    for ref_fasta_file in ref_fasta_filelist:
        ref_fasta_path = os.path.join(kallisto_ref_dir, ref_fasta_file)
        genome_id = ref_fasta_file.split(".fna")[0]
        index_file = genome_id + ".idx"
        index_path = os.path.join(kallisto_index_dir, index_file)
        kallisto_index_args = ["kallisto", "index", "-i", index_path, ref_fasta_path]
        subprocess.run(kallisto_index_args)
    return


def run_kallisto_quant(Hawkey2022_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, results_dir):
    index_list = [x for x in os.listdir(kallisto_index_dir) if x.endswith(".idx")]
    for index_file in index_list:
        index_path = os.path.join(kallisto_index_dir, index_file)
        genome_id = index_file.split(".idx")[0]
        SRA_id = Hawkey2022_genomeID_to_SRA_ID_dict[genome_id]
        if SRA_id == "NA":
            continue
        else:
            read_file1 = SRA_id + "_1.fastq"
            read_file2 = SRA_id + "_2.fastq"
            read_path1 = os.path.join(SRA_data_dir, read_file1)
            read_path2 = os.path.join(SRA_data_dir, read_file2)
            output_path = os.path.join(results_dir, genome_id)
            ## run with 10 threads by default.
            kallisto_quant_args = ["kallisto", "quant", "-t", "10", "-i", index_path, "-o", output_path, "-b", "100", read_path1, read_path2]
            subprocess.run(kallisto_quant_args)
    return


def make_genome_to_SRA_dict(Hawkey2022_metadata_csv):
    genome_to_SRA_dict = dict()
    with open(Hawkey2022_metadata_csv, "r") as csv_fh:
        for i, line in enumerate(csv_fh):
            if i == 0: continue ## skip the header.
            line = line.strip()
            NCBI_bioproject, ReferenceGenome, BioSample, SRA_Data = line.split(',')
            genome_id = ReferenceGenome.split(".gbff.gz")[0]
            genome_to_SRA_dict[genome_id] = SRA_Data
    return genome_to_SRA_dict


def parse_metadata_in_header(target_id):
    fields = target_id.split("|")
    SeqID = fields[0].split("=")[-1]
    SeqType = fields[1].split("=")[-1]
    locus_tag = fields[2].split("=")[-1]
    ## convert underscores back into spaces.
    product = fields[3].split("=")[-1].replace("_", " ")
    metadata_tuple = (SeqID, SeqType, locus_tag, product)
    return(metadata_tuple)


def estimate_chr_plasmid_copy_numbers(genecount_tsv_path):
    genome_dict = dict()
    ## keys are SeqIDs.
    ## values are a dict: {SeqType: "chromosome", total_length: 10000, total_est_counts: 100}
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_metadata_in_header(target_id)
            if SeqID in genome_dict:
                genome_dict[SeqID]["total_length"] += float(length)
                genome_dict[SeqID]["total_est_counts"] += float(est_counts)
            else: ## Initialize the dictionary.
                genome_dict[SeqID] = {"SeqType" : SeqType, "total_length" : float(length), "total_est_counts": float(est_counts)}
    coverage_dict = dict()
    ##keys are seq_ids, value is (SeqType, coverage) pair.
    chromosome_coverage = -1
    for SeqID, replicon_dict in genome_dict.items():
        coverage = replicon_dict["total_est_counts"]/replicon_dict["total_length"]
        coverage_dict[SeqID] = (replicon_dict["SeqType"], coverage)
        if replicon_dict["SeqType"] == "chromosome":
            chromosome_coverage = coverage
    ## now normalize by chromosome coverage to get copy number estimates.
    copy_number_dict = dict()
    for SeqID, value_tuple in coverage_dict.items():
        seqtype, coverage = value_tuple
        copy_number_dict[SeqID] = (seqtype, coverage/chromosome_coverage)
    return(copy_number_dict)


def measure_Hawkey2022_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    AnnotationAccession, SeqID, SeqType, CopyNumber
    """
    AnnotationAccessionVec = []
    SeqIDVec = []
    SeqTypeVec = []
    CopyNumberVec = []
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        ## I probably should have trimmed the '_genomic' suffix in an earlier step.
        annotation_accession = genomedir.split("_genomic")[0]
        genome_quantfile_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
        copy_number_dict = estimate_chr_plasmid_copy_numbers(genome_quantfile_path)
        for SeqID, value_tuple in copy_number_dict.items():
            seqtype, coverage = value_tuple
            AnnotationAccessionVec.append(annotation_accession)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            CopyNumberVec.append(coverage)

    assert len(AnnotationAccessionVec) == len(SeqIDVec) == len(SeqTypeVec) == len(CopyNumberVec)
    ## now write the copy number data to file.
    with open(copy_number_csv_file, "w") as outfh:
        header = "AnnotationAccession,SeqID,SeqType,CopyNumber"
        outfh.write(header + "\n")
        for i in range(len(AnnotationAccessionVec)):
            outfh.write(AnnotationAccessionVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + str(CopyNumberVec[i]) + "\n")
    return


################################################################################################
def isARG(product_annotation):
    chloramphenicol_keywords = "chloramphenicol|Chloramphenicol"
    tetracycline_keywords = "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
    MLS_keywords = "macrolide|lincosamide|streptogramin"
    multidrug_keywords = "Multidrug resistance|multidrug resistance|antibiotic resistance"
    beta_lactam_keywords = "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\S*"
    glycopeptide_keywords = "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
    polypeptide_keywords = "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
    diaminopyrimidine_keywords = "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
    sulfonamide_keywords = "sulfonamide|Sul1|sul1|sulphonamide"
    quinolone_keywords = "quinolone|Quinolone|oxacin|qnr|Qnr"
    aminoglycoside_keywords = "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
    macrolide_keywords = "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
    antimicrobial_keywords = "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\S*"
    ARG_regex = "|".join([chloramphenicol_keywords, tetracycline_keywords,
                          MLS_keywords, multidrug_keywords, beta_lactam_keywords,
                          glycopeptide_keywords, polypeptide_keywords, diaminopyrimidine_keywords,
                          sulfonamide_keywords, quinolone_keywords, aminoglycoside_keywords,
                          macrolide_keywords, antimicrobial_keywords])
    if re.search(ARG_regex, product_annotation): return True
    return False


def estimate_ARG_copy_numbers(genecount_tsv_path):

    chromosomal_gene_length = 0.0
    chromosomal_gene_est_counts = 0.0

    ARG_coverage_dict = dict()
    ## get the chromosomal gene coverage, and get the coverage for all ARGs.
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_metadata_in_header(target_id)
            if SeqType == "chromosome":
                chromosomal_gene_length += float(length)
                chromosomal_gene_est_counts += float(est_counts)
            if isARG(product):
                coverage = float(est_counts) / float(length)
                ARG_coverage_dict[locus_tag] = (SeqID, SeqType, product, coverage)
    ## NOTE: GCF_026154285.1_ASM2615428v1 did not have any reads pseudoalign.
    ## Return an empty dict() when nothing aligns to the chromosome.
    if chromosomal_gene_length == 0:
        print("WARNING: no reads pseudoaligned in file: ", genecount_tsv_path)
        print("estimate_ARG_copy_numbers is returning an empty dict.")
        return(dict())
    chromosome_coverage = chromosomal_gene_est_counts/chromosomal_gene_length
    ## now normalize by chromosome coverage to get copy number estimates.
    ARG_copy_number_dict = dict()
    for locus_tag, value_tuple in ARG_coverage_dict.items():
        my_SeqID, my_SeqType, my_product, my_coverage = value_tuple
        ARG_copy_number_dict[locus_tag] = (my_SeqID, my_SeqType, my_product, my_coverage/chromosome_coverage)
    return(ARG_copy_number_dict)


def measure_Hawkey2022_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    AnnotationAccession, SeqID, SeqType, locus_tag, product, CopyNumber
    """
    AnnotationAccessionVec = []
    SeqIDVec = [] ## this is for the replicon.
    SeqTypeVec = []
    LocusTagVec = []
    ProductVec = []
    CopyNumberVec = []
    
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        ## I probably should have trimmed the '_genomic' suffix in an earlier step.
        annotation_accession = genomedir.split("_genomic")[0]
        genome_quantfile_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
        ARG_copy_number_dict = estimate_ARG_copy_numbers(genome_quantfile_path)
        for locus_tag, value_tuple in ARG_copy_number_dict.items():
            SeqID, seqtype, product, copy_number = value_tuple
            AnnotationAccessionVec.append(annotation_accession)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            LocusTagVec.append(locus_tag)
            ProductVec.append(product)
            CopyNumberVec.append(copy_number)

    assert len(AnnotationAccessionVec) == len(SeqIDVec) == len(SeqTypeVec) == len(LocusTagVec) == len(ProductVec) == len(CopyNumberVec)
    ## now write the ARG copy number data to file.
    with open(ARG_copy_number_csv_file, "w") as outfh:
        header = "AnnotationAccession,SeqID,SeqType,locus_tag,product,CopyNumber"
        outfh.write(header + "\n")
        for i in range(len(AnnotationAccessionVec)):
            outfh.write(AnnotationAccessionVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + LocusTagVec[i] + "," + ProductVec[i] + "," + str(CopyNumberVec[i]) + "\n")
    return


def main():

    SRA_data_dir = "../data/Hawkey2022-SRA-data/"
    SRA_accession_list_file = os.path.join(SRA_data_dir, "Hawkey2022-SRA-accessions.txt")
    refgenomes_dir = "../data/Hawkey2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA646837/"
    kallisto_ref_dir = "../results/Hawkey2022_kallisto_references/"
    kallisto_index_dir = "../results/Hawkey2022_kallisto_indices/"
    kallisto_quant_results_dir = "../results/Hawkey2022_kallisto_quantification"
    Hawkey2022_metadata_csv = os.path.join(SRA_data_dir, "Hawkey2022-genome-metadata.csv")
    Hawkey2022_genomeID_to_SRA_ID_dict = make_genome_to_SRA_dict(Hawkey2022_metadata_csv)
    copy_number_csv_file = "../results/Hawkey2022_chromosome_plasmid_copy_numbers.csv"
    ARG_copy_number_csv_file = "../results/Hawkey2022_ARG_copy_numbers.csv"
    
    download_fastq_reads(SRA_data_dir, SRA_accession_list_file)
    make_Hawkey_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_dir)
    make_Hawkey_kallisto_indices(kallisto_ref_dir, kallisto_index_dir)
    run_kallisto_quant(Hawkey2022_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, kallisto_quant_results_dir)
    measure_Hawkey2022_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file)
    measure_Hawkey2022_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file)
    return


main()
