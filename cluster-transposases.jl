"""
cluster-transposases.jl by Rohan Maddamsetti.

This script merges transposases associated with duplicated ARGs in E. coli that are 99% similar
by comparison to a reference transposase sequence (the most common one in the cluster).

This script also merges transposases that are highly similar. for plotting their distribution across taxa.

The input for this script is produced by ARG-duplication-analysis.R.
For now, the output for this script is also used by ARG-duplication-analysis.R to make figures.

usage: julia cluster-transposases.jl
"""

using DataFrames, DataFramesMeta, CSV, BioSequences, BioAlignments, FASTX


mutable struct EcoliTransposonCluster
    ## the most common sequence is set as the reference.
    ref::LongAA
    product_annotation::String
    counts::Int64
end


function ClusterEcoliTransposases(transposase_df)
    ## create an affine gap scoring model
    affinegap = AffineGapScoreModel(
        match=1,
        mismatch=-1,
        gap_open=-1,
        gap_extend=-1
    )
    
    clusters = Vector{EcoliTransposonCluster}()
    ## critical assumption: the sequences are read in sorted order
    ## by the number of transposase counts (most to least)
    for (seq, product_annot, count) in zip(transposase_df.sequence, transposase_df.product, transposase_df.count)
        prot_seq = LongAA(seq)
        ## scale the number of mismatches by the length of the sequence (99% similarity)
        max_mismatches = floor(length(prot_seq) * 0.01)
        matched_cluster = false
        for cluster in clusters
            res = pairalign(GlobalAlignment(), prot_seq, cluster.ref, affinegap)
            aln = alignment(res)
            mismatches = count_mismatches(aln)
            insertions = count_insertions(aln)
            deletions = count_deletions(aln)
            if (mismatches + insertions + deletions) <= max_mismatches
                matched_cluster = true
                cluster.counts += count
                break
            end
        end
        if matched_cluster == false
            newCluster = EcoliTransposonCluster(prot_seq, product_annot, count)
            push!(clusters, newCluster)
        end
    end
    return clusters
end


function WriteEcoliClustersToFile(clusters, outfile)
    outfh = open(outfile, "w")
    header = "product,sequence,count"
    println(outfh, header)
    for cluster in clusters
        line = string(cluster.product_annotation) * "," * string(cluster.ref) * "," * string(cluster.counts)
        println(outfh, line)
    end
    close(outfh)
end


mutable struct GenusTransposonCluster
    ## the most common sequence is set as the reference.
    ref::LongAA
    product_annotation::String
    genus_to_counts::Dict{String,Int64}
end


function ClusterGenusTransposases(all_transposase_df)
    
    ## we need to sort the transposases by how often they are observed in the data,
    ## so that we use the most frequent sequence as the reference for a cluster.
    sorted_df = @chain all_transposase_df begin
        groupby([:sequence, :Genus])
        @combine begin
            :count = length(:sequence)
        end
        @orderby -:count ## put in descending order.
    end

    ## make a product_annotation column for sorted_transposase_df
    ## by greedily searching all_transposase_df for matches.
    sorted_product_annotation_vec = []
    for seq in sorted_df.sequence
        for (matching_seq, product_annot) in zip(all_transposase_df.sequence, all_transposase_df.product_annotation)
            if (seq == matching_seq)
                sorted_product_annotation_vec = vcat(sorted_product_annotation_vec, product_annot)
                break
            else
                continue
            end
        end
    end
    @assert length(sorted_product_annotation_vec) == length(sorted_df.sequence)
    ## only add the new product_annotation column.
    sorted_df.product_annotation = sorted_product_annotation_vec
    
    ## create an affine gap scoring model
    affinegap = AffineGapScoreModel(
        match=1,
        mismatch=-1,
        gap_open=-1,
        gap_extend=-1
    )

    clusters = Vector{GenusTransposonCluster}()
    ## critical assumption: the sequences are read in sorted order
    ## by the number of transposase counts (most to least)
    for (seq, genus, product_annot, count) in zip(sorted_df.sequence, sorted_df.Genus, sorted_df.product_annotation, sorted_df.count)
        prot_seq = LongAA(seq)
        ## scale the number of mismatches by the length of the sequence (99% similarity)
        max_mismatches = floor(length(prot_seq) * 0.01)
        matched_cluster = false
        for cluster in clusters
            res = pairalign(GlobalAlignment(), prot_seq, cluster.ref, affinegap)
            aln = alignment(res)
            mismatches = count_mismatches(aln)
            insertions = count_insertions(aln)
            deletions = count_deletions(aln)
            if (mismatches + insertions + deletions) <= max_mismatches
                matched_cluster = true
                if haskey(cluster.genus_to_counts, genus)
                    cluster.genus_to_counts[genus] += count
                else
                    cluster.genus_to_counts[genus] = count
                end
                break
            end
        end
        if matched_cluster == false
            new_genus_to_counts = Dict{String,Int64}([(genus, count)])
            newCluster = GenusTransposonCluster(prot_seq, product_annot, new_genus_to_counts)
            push!(clusters, newCluster)
        end
    end
    return clusters
end


function WriteGenusClustersToFile(genus_clusters, outfile)
    outfh = open(outfile, "w")
    header = "product, Genus, sequence,count"
    println(outfh, header)
    for cluster in genus_clusters
        for (genus, counts) in cluster.genus_to_counts
        line = string(cluster.product_annotation) * "," * string(genus) * "," * string(cluster.ref) * "," * string(counts)
        println(outfh, line)
        end
    end
    close(outfh)
end


function main()
        
    ## let's start off by clustering the E. coli transposons associated with duplicated ARGs.
    Ecoli_transposase_path = "../results/Ecoli-transposases-in-dup-regions-with-ARGs.csv"
    Ecoli_transposase_df = CSV.read(Ecoli_transposase_path, DataFrame)
    
    ## cluster the transposases.
    ecoli_clusters = ClusterEcoliTransposases(Ecoli_transposase_df)
    ## and write to file
    Ecoli_cluster_outfile = "../results/merged_Ecoli-transposases-in-dup-regions-with-ARGs.csv"
    WriteEcoliClustersToFile(ecoli_clusters, Ecoli_cluster_outfile)

    ## there are 685 E. coli transposases associated with duplicated ARGs.
    total_counts = 0
    for x in ecoli_clusters
        total_counts += x.counts
    end
    println("E. coli transposases associated with duplicated ARGs:  " * string(total_counts))
    ## IS26 (at 99% similarity threshold) accounts for 418 of these.
    println("The top transposase (IS26) is associated with duplicated ARGs this many times: " * string(ecoli_clusters[1].counts))

    ## Now, let's cluster all transposases associated with duplicated ARGs, using a 99% similarity threshold.
    all_transposase_path = "../results/transposases-in-dup-regions-with-ARGs.csv"
    all_transposase_df = CSV.read(all_transposase_path, DataFrame)
    ## cluster the transposases.
    genus_clusters = ClusterGenusTransposases(all_transposase_df)
    ## and write to file
    genus_cluster_outfile = "../results/merged_transposases-in-dup-regions-with-ARGs.csv"
    WriteGenusClustersToFile(genus_clusters, genus_cluster_outfile)

end


main()
