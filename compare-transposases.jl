"""
compare-transposases.jl by Rohan Maddamsetti.

Align the transposons in Yi's library with the transposases associated with
 duplicated ARGs, to see if we have any of these on hand in the lab.

All comparisons are done on amino acid sequences.

CRITICAL NOTE:
I used the ISFinder database rather than using the sequences listed on Benchling.
So, there is the potential for inconsistencies between the actual sequences in
Yi's transposon library, and the sequence matching the given ID in ISFinder.

I found a couple such cases already. But for my current intents and purposes,
this doesn't really matter so much.


Input files:
=================

-- A table of the native transposases in Yi's library:
../data/YiYao-native-transposon-library-IS-Finder-seqs.csv

-- A table (generated by ARG-duplication-analysis.R) of transposases in
regions of duplicated genes that contain ARGs in the complete genomes on Genbank:
../results/transposases-in-dup-regions-with-ARGs.csv


Usage examples:
==================
julia compare-transposases.jl

"""

using DataFrames, DataFramesMeta, CSV, BioSequences, BioAlignments, FASTX, CodecZlib


function ClusterTransposases(yao_library_dict, duplicated_transposase_df, max_mismatches, affinegap)

    ## let's keep things really simple at first.    
    clusters = Dict((k, Vector()) for k ∈ keys(yao_library_dict))

    for seq in duplicated_transposase_df.sequence
        AAseq = LongAA(seq)
        for yaoID in keys(yao_library_dict)
            refseq = LongAA(yao_library_dict[yaoID])
            res = pairalign(GlobalAlignment(), AAseq, refseq, affinegap)
            aln = alignment(res)
            mismatches = count_mismatches(aln)
            insertions = count_insertions(aln)
            deletions = count_deletions(aln)
            if (mismatches + insertions + deletions) <= max_mismatches
                push!(clusters[yaoID], AAseq)
                break
            end
        end
    end
    return clusters
end


function main()

    max_mismatches = 3 ## perfect matches required.

    """ The table of the native transposases in Yi's library. """
    YaoLibraryCSV = "../data/YiYao-native-transposon-library-IS-Finder-seqs.csv"
    yao_library_df = DataFrame(CSV.File(YaoLibraryCSV))

    """" The table (generated by ARG-duplication-analysis.R) of transposases in
    regions of duplicated genes that contain ARGs in the complete genomes on Genbank. """

    DuplicatedTransposasesCSV = "../results/transposases-in-dup-regions-with-ARGs.csv"
    duplicated_transposase_df = DataFrame(CSV.File(DuplicatedTransposasesCSV))
    
    ## create an affine gap scoring model
    affinegap = AffineGapScoreModel(
        match=1,
        mismatch=-1,
        gap_open=-1,
        gap_extend=-1
    )

    ## create a lookup table between Yao Catalogy ID and AA sequences.
    yao_library_dict = Dict((k, v) for (k, v) ∈ zip(yao_library_df.Catalog_number, yao_library_df.sequence))
    
    ## cluster the transposases.
    ## This takes 106 seconds.
    @time clusters = ClusterTransposases(yao_library_dict, duplicated_transposase_df, max_mismatches, affinegap)
    """"
     perfect matches to the following transposases in Yi's library:
     pB112 (IS6-IS26). ISFinder reports that this element is associated with many different ARGs in nature.
    Also see this mBio paper: "Movement of IS26-Associated Antibiotic Resistance Genes Occurs
    via a Translocatable Unit That Includes a Single IS26 and Preferentially Inserts Adjacent to Another IS26"

     pB123 (IS256-ISEc58).
     pB145 (ISNCY-IS1202-ISKpn21).
     pB146 (Tn3-Tn2). This one has a passenger ARG in ISFinder.
    
     pB111 (IS5-IS903) has a 3-mismatch neighbor among the duplicated transposases.
    """
    
    ## get matching sequences among the duplicated transposases.
    ## union gets rid of duplicate elements in the array.
    matched_sequences = union(ExtractMatchedSequencesFromTransposaseClusters(clusters))
    ## now subset duplicated_transposases_df based on the matches to get the metadata.
    matched_duplicated_transposase_df = @rsubset(duplicated_transposase_df, :sequence in matched_sequences)

    outfile = "../results/duplicated-transposases-matching-Yao-library-transposons.csv"
    ## if write matches to file.
    CSV.write(outfile, matched_duplicated_transposase_df)
end


main()