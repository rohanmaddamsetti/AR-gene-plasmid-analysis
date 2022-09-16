from unittest import TestCase, main
from io import StringIO

from genomediff import Metadata, GenomeDiff
from genomediff.parser import GenomeDiffParser
from genomediff.records import Record


class ParserTestCase(TestCase):
    def test_parse(self):
        file = StringIO("""
#=GENOME_DIFF	1.0
#=AUTHOR test
SNP	1	23423	NC_000913	223	A	gene_name=mhpE
RA	2		NC_000913	223	0	G	A	frequency=0.1366
                        """.strip())
        p = GenomeDiffParser(fsock=file)
        self.assertEqual([
                             Metadata('GENOME_DIFF', '1.0'),
                             Metadata('AUTHOR', 'test'),
                             Record('SNP', 1, parent_ids=[23423], new_seq='A', seq_id='NC_000913', position=223, gene_name='mhpE'),
                             Record('RA', 2, new_base='A', frequency=0.1366, position=223, seq_id='NC_000913',
                                    insert_position=0,
                                    ref_base='G')],
                         list(p)
        )

    def test_parse_dot_missing_parent_ids(self):
        file = StringIO("""
#=GENOME_DIFF	1.0
#=AUTHOR test
SNP	1	23423	NC_000913	223	A	gene_name=mhpE
RA	2	.	NC_000913	223	0	G	A	frequency=0.1366
                        """.strip())
        p = GenomeDiffParser(fsock=file)
        self.assertEqual([
                             Metadata('GENOME_DIFF', '1.0'),
                             Metadata('AUTHOR', 'test'),
                             Record('SNP', 1, parent_ids=[23423], new_seq='A', seq_id='NC_000913', position=223, gene_name='mhpE'),
                             Record('RA', 2, new_base='A', frequency=0.1366, position=223, seq_id='NC_000913',
                                    insert_position=0,
                                    ref_base='G')],
                         list(p)
        )


class GenomeDiffTestCase(TestCase):
    def test_document(self):
        file = StringIO("""
#=GENOME_DIFF	1.0
#=AUTHOR test
SNP	1	23423	NC_000913	223	A
RA	2		NC_000913	223	0	G	A
                        """.strip())

        document = GenomeDiff.read(file)

        self.assertEqual({'AUTHOR': 'test', 'GENOME_DIFF': '1.0'}, document.metadata)

        snp_record = Record('SNP', 1, document, [23423], seq_id='NC_000913', new_seq='A', position=223)
        ra_record = Record('RA', 2, document, None, position=223, seq_id='NC_000913', insert_position=0, new_base='A',
                           ref_base='G')

        self.assertEqual([snp_record], document.mutations)
        self.assertEqual([ra_record], document.evidence)
        self.assertEqual(snp_record, document[1])
        self.assertEqual(ra_record, document[2])


class RecordTestCase(TestCase):
    def test_simple(self):
        snp_record = Record('SNP', 1, parent_ids=[23423], seq_id='NC_000913', new_seq='A', position=223, test='more')

        self.assertEqual('SNP', snp_record.type)
        self.assertEqual(1, snp_record.id)
        self.assertEqual('A', snp_record.new_seq)
        self.assertEqual('more', snp_record.test)


class ParentResolveTestCase(TestCase):
    def test_resolve(self):
        file = StringIO("""
#=GENOME_DIFF	1.0
#=AUTHOR test
SNP	1	2	NC_000913	223	A
RA	2		NC_000913	223	0	G	A
                        """.strip())
        document = GenomeDiff.read(file)
        self.assertEqual(document[1].parents, [document[2]])

class RecordComparisonTestCase(TestCase):
    def test_cmp1(self):
        file1 = StringIO("""
#=GENOME_DIFF	1.0
#=CREATED	20:02:17 23 Jan 2019
#=PROGRAM	breseq 0.33.2 
#=COMMAND	breseq -r LCA.gff3 sequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R1.fastq.gz sequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R2.fastq.gz sequence-data/ZDBp889_reads.fastq -o consensus/ZDBp889
#=REFSEQ	LCA.gff3
#=READSEQ	sequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R1.fastq.gz
#=READSEQ	sequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R2.fastq.gz
#=READSEQ	sequence-data/ZDBp889_reads.fastq
#=CONVERTED-BASES	644779377
#=CONVERTED-READS	14448149
#=INPUT-BASES	645034321
#=INPUT-READS	14455411
#=MAPPED-BASES	602854657
#=MAPPED-READS	13788351
SNP	1	34	REL606	72313	C
        """.strip())

        document1 = GenomeDiff.read(file1)
        
        file2 = StringIO("""
#=GENOME_DIFF	1.0
#=CREATED	16:49:49 23 Jan 2019
#=PROGRAM	breseq 0.33.2 
#=COMMAND	breseq -r LCA.gff3 sequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R1.fastq.gz sequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R2.fastq.gz -o consensus/ZDB67
#=REFSEQ	LCA.gff3
#=READSEQ	sequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R1.fastq.gz
#=READSEQ	sequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R2.fastq.gz
#=CONVERTED-BASES	114566968
#=CONVERTED-READS	419781
#=INPUT-BASES	114567554
#=INPUT-READS	419783
#=MAPPED-BASES	92472620
#=MAPPED-READS	339813
SNP	1	12	REL606	72313	C
        """.strip())

        document2 = GenomeDiff.read(file2)
        self.assertEqual(document1.mutations,document2.mutations)


    def test_cmp2(self):
        file1 = StringIO("""
#=GENOME_DIFF	1.0
SNP	1	12	REL606	72313	C	aa_new_seq=G	aa_position=92	aa_ref_seq=D	codon_new_seq=GGC	codon_number=92	codon_position=2	codon_ref_seq=GAC	gene_name=araA	gene_position=275	gene_product=L-arabinose isomerase	gene_strand=<	genes_overlapping=araA	locus_tag=ECB_00064	locus_tags_overlapping=ECB_00064	mutation_category=snp_nonsynonymous	position_end=72313	position_start=72313	snp_type=nonsynonymous	transl_table=11
        """.strip())

        document1 = GenomeDiff.read(file1)
        
        file2 = StringIO("""
#=GENOME_DIFF	1.0
SNP	1	34	REL606	72313	C	aa_new_seq=G	aa_position=92	aa_ref_seq=D	codon_new_seq=GGC	codon_number=92	codon_position=2	codon_ref_seq=GAC	gene_name=araA	gene_position=275	gene_product=L-arabinose isomerase	gene_strand=<	genes_overlapping=araA	locus_tag=ECB_00064	locus_tags_overlapping=ECB_00064	mutation_category=snp_nonsynonymous	position_end=72313	position_start=72313	snp_type=nonsynonymous	transl_table=11
        """.strip())

        document2 = GenomeDiff.read(file2)
        self.assertEqual(document1.mutations,document2.mutations)

class RecordSatisfiesTestCase(TestCase):
    def test_satisfies1(self):

        file1 = StringIO("""
#=GENOME_DIFF	1.0
SNP	1	12	REL606	72313	C	aa_new_seq=G	aa_position=92	aa_ref_seq=D	codon_new_seq=GGC	codon_number=92	codon_position=2	codon_ref_seq=GAC	gene_name=araA	gene_position=275	gene_product=L-arabinose isomerase	gene_strand=<	genes_overlapping=araA	locus_tag=ECB_00064	locus_tags_overlapping=ECB_00064	mutation_category=snp_nonsynonymous	position_end=72313	position_start=72313	snp_type=nonsynonymous	transl_table=11
        """.strip())

        document1 = GenomeDiff.read(file1)
        
        file2 = StringIO("""
#=GENOME_DIFF	1.0
SNP	1	34	REL606	72313	C	aa_new_seq=G	aa_position=92	aa_ref_seq=D	codon_new_seq=GGC	codon_number=92	codon_position=2	codon_ref_seq=GAC	gene_name=araA	gene_position=275	gene_product=L-arabinose isomerase	gene_strand=<	genes_overlapping=araA	locus_tag=ECB_00064	locus_tags_overlapping=ECB_00064	mutation_category=snp_nonsynonymous	position_end=72313	position_start=72313	snp_type=nonsynonymous	transl_table=11
        """.strip())

        document2 = GenomeDiff.read(file2)
        
        self.assertEqual(document1.mutations[0].satisfies("gene_name==araA"),True)
        self.assertEqual(document1.mutations[0].satisfies("mutation_category==snp_nonsynonymous"),True)
        self.assertEqual(document1.mutations[0].satisfies("position_start>=50.5"),True)
        self.assertEqual(document1.mutations[0].satisfies("position_start<=10000000.2"),True)

        document2.remove("gene_name==rrlA")
        print()
        print(document1.mutations)
        print()
        print(document2.mutations)
        self.assertEqual(document1.mutations,document2.mutations)
        document2.remove("gene_name==araA")
        self.assertEqual(document2.mutations,[])
        
if __name__ == '__main__':
    main()
