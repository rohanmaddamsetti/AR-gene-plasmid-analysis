### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1ba467c8-5a15-4342-adca-486c970a0c0f
using CSV, DataFrames, Plots, LightGraphs, SparseArrays, LinearAlgebra, Distances, Clustering, MinHash

# ╔═╡ 73979b34-caea-11eb-3431-8d722c6416b9
md""" ### HGT-network-notebook.jl by Rohan Maddamsetti

## Fast gene flow network inference from copy number variation in bacterial genomes

#### Rohan Maddamsetti  

The initial goal of this notebook is to draft matrix and graph visualizations of how identical protein sequences are shared across strains in my dataset.

I will do analyses that focus just on duplicated genes, as well as singleton genes. I I will then repeat these analyses, focusing on genes associated with mobile genetic elements.

"""

# ╔═╡ 0c0a7d17-44a2-4ec7-aa83-b1d52920cd2f
md""" IMPORTANT NOTE!!! Use MinHash algorithm (Mash or Sourmash software) in order to rapidly estimate distances and cluster genes and strains. I need something that scales!!! Since I'm not aware of any mature implementations of MinHash in Julia, I rewrote this code in python, using the DataSketch module to use MinHash. """

# ╔═╡ cf2d881b-0502-4ee0-8b49-94c33427c263
md""" let's look at the duplicate proteins, first."""

# ╔═╡ ec356b3d-18a3-4977-8fba-2ca9ee7524a7
begin
	
	dup_protein_fpath = "../results/AR-gene-duplication/duplicate-proteins.csv"
	dup_protein_file = CSV.File(dup_protein_fpath)
	dup_protein_df = DataFrame(dup_protein_file)
	
end

# ╔═╡ 2afb1c61-50a2-4521-8018-9cc5e9b014b6
md""" 
now, let's visualize a graph between these duplicate sequences, and the genomes that they come from. We can use a matrix "heatmap" to do this visualization, as well as a network visualization of the graph.
"""

# ╔═╡ 25f6e9b8-a433-4eb1-92f2-91aad7c65e87
function ProteinToAccession(protein_df)
	## make a dictionary of protein sequence to Annotation_Accession.
	protseq_to_accession = Dict()
	
	for row in eachrow(protein_df)
		accession = row.Annotation_Accession
		sequence = row.sequence

		if haskey(protseq_to_accession, sequence)
			push!(protseq_to_accession[sequence], accession)
		else
			protseq_to_accession[sequence] = [accession]
		end
	end
	
	return protseq_to_accession
end

# ╔═╡ d492149c-2fe7-4992-a0bb-3c81f27b1b30
function SharedProteinToAccession(protseq_to_accession)
	## filter for sequences that are present in multiple genomes.
	HGT_protseq_to_accession = Dict()
	for (k,v) in protseq_to_accession
		if length(v) > 1
			HGT_protseq_to_accession[k] = v
		end
	end
	return HGT_protseq_to_accession
end

# ╔═╡ 463fea4c-ccd1-4b70-9e38-0d8f31a523e1
function AccessionToProtein(protein_df)
	## make a dictionary of Annotation_Accession to protein sequence.
	accession_to_protseq = Dict()
	
	for row in eachrow(protein_df)
		accession = row.Annotation_Accession
		sequence = row.sequence

		if haskey(accession_to_protseq, accession)
			push!(accession_to_protseq[accession], sequence)
		else
			accession_to_protseq[accession] = [sequence]
		end
	end
	
	return accession_to_protseq
end

# ╔═╡ f96f9c72-f05a-4998-8588-555120dc61e9
function FilterAccessionToProteinDictForSharedSequences(accession_to_protseq, HGT_protseq_to_accession)
	## filter for sequences that are present in multiple genomes.
	filtered_accession_to_protseq = Dict()
	
	for (accession,protvec) in accession_to_protseq
		for prot in protvec
			if prot in keys(HGT_protseq_to_accession)
				filtered_accession_to_protseq[accession] = protvec
			end
		end
	end

	return filtered_accession_to_protseq
end

# ╔═╡ b1509dee-b5a0-4dd2-a1eb-60a82bf21866
## make a dictionary of protein sequence to Annotation_Accession.
HGT_dup_protseq_to_accession = dup_protein_df |>
ProteinToAccession |> SharedProteinToAccession

# ╔═╡ 6b5ef553-0bb4-47d1-b6cf-c6fb898675f3
md"""

A good heuristic for relatedness is whether two organisms have identical copies of the tuf gene encoding the Tu elongation factor. If two organisms do not share identical tuf sequences, but they have other identical proteins, then it is likely that those identical proteins stem from HGT rather than vertical descent. The same heuristic can be extended to conserved phylogenetic marker genes, etc.

"""

# ╔═╡ 86fc624f-4874-4aa4-acd7-6d70ff9c927f
HGT_dup_protein_df = filter(row -> haskey(HGT_dup_protseq_to_accession,row.sequence), dup_protein_df)

# ╔═╡ d63763c7-ce0d-4f7d-b7d9-2e9b5a8a5fa9
## make a dictionary of Annotation_Accession to protein sequence.
HGT_dup_accession_to_protseq = AccessionToProtein(HGT_dup_protein_df)

# ╔═╡ e75548f9-60ba-41c5-ad83-c3a1061b4177
plasmid_HGT_dup_protein_df = filter(row -> row.plasmid_count > 1, HGT_dup_protein_df)

# ╔═╡ 8cfe264a-187f-4c63-b5e0-80cdaa284896
md"""

let's actually make some graphs now.


IDEA: in my plots, I can place "Gene" nodes based on their genetic distance to each other. In addition, "Strain" nodes can also have a genetic distance as well. How can I make use of this information as well?

first, make a matrix plot, where strains are rows, and genes are columns.

then, plot a bipartite graph, where the two node types are strains and genes.

"""

# ╔═╡ d495cb9f-02e0-45eb-b5d8-0a091faf9d8e
function ProteinToAccessionMatrix(protseq_to_accession_dict)
	
	## make a matrix as the data structure for the plot.
	## rows are the protein sequences.
	## columns are the Accessions.
	
	accession_set = Set()
	seq_set = Set()
	
	for (seq, accession_list) in protseq_to_accession_dict
		push!(seq_set, seq)
		for a in accession_list
			push!(accession_set, a)
		end
	end
	
	n_accessions = length(accession_set)
	n_seqs = length(seq_set)
	
	sorted_accessions = sort([x for x in accession_set])
	sorted_seqs = sort([x for x in seq_set])
	
	accession_to_index = Dict(x => i for (i,x) in enumerate(sorted_accessions))
	seq_to_index = Dict(x => i for (i,x) in enumerate(sorted_seqs))
	
	accessions_seq_matrix = zeros((n_accessions, n_seqs))
	
	for (seq, accession_list) in protseq_to_accession_dict
		seq_idx = seq_to_index[seq]
		for a in accession_list
			accession_idx = accession_to_index[a]
			accessions_seq_matrix[accession_idx, seq_idx] = 1
		end
	end
	
	return accessions_seq_matrix
	
end

# ╔═╡ 196b9ad6-5251-4534-bb7b-185130589ffb
md"""  Cluster the rows, and cluster the columns.
 sort and update the seq to index and accession to index mappings
 by getting the permutation of columns and rows returned by the clustering.

 Speed up these pairwise calculations by using MinHash! See implementation at:
https://jakobnissen.github.io/MinHash.jl"""

# ╔═╡ 2c068fa2-83a9-466a-be31-3d480e18aa21
function unionlength(x::MinHashSketch, y::MinHashSketch)

    xi, yi = 1, 1
    n = 0

    @inbounds while (xi ≤ length(x.hashes)) & (yi ≤ length(y.hashes))
        xv, yv = x.hashes[xi], y.hashes[yi]
        n += 1
        xi += xv ≤ yv
        yi += yv ≤ xv
    end

    return n

end

# ╔═╡ 3736e126-ebb4-4949-83a9-14d7922dc3d3
function unionlength_vec(hash_vec)
	union_count_matrix = zeros(Int32, length(hash_vec), length(hash_vec))
	for i in 1:length(hash_vec)
		for j in 1:length(hash_vec)
			if i >= j
				union_count_matrix[i,j] = unionlength(hash_vec[i], hash_vec[j])
			end
		end
	end
	
	return union_count_matrix
			
end

# ╔═╡ c8c5417f-e7c4-4e68-874f-5ed9d15139da
function SortAccessionSeqMatrix(protseq_to_accession, accession_to_protseq)
	
	accessions_seq_matrix = ProteinToAccessionMatrix(protseq_to_accession)
	
	## use MinHash to calculate the Jaccard similarity of the columns (genes).
	genome_hashes = [MinHash.sketch(accession_to_protseq[k], 500) for k in sort([x for x in keys(accession_to_protseq)])]
	
	genome_intersection_matrix = intersectionlength(genome_hashes)
	
	## use MinHash to calculate the Jaccard similarity of the columns (genes).
	gene_hashes = [MinHash.sketch(protseq_to_accession[k], 500) for k in sort([x for x in keys(protseq_to_accession)])]
	
	gene_intersection_matrix = intersectionlength(gene_hashes)
	
	genome_union_matrix = unionlength_vec(genome_hashes)
	gene_union_matrix = unionlength_vec(gene_hashes)
	
	genome_distance_matrix = 1 .- (genome_intersection_matrix./genome_union_matrix) - I ## subtract identity matrix
	replace!(genome_distance_matrix, NaN=>0)
	
	gene_distance_matrix = 1 .- (gene_intersection_matrix./gene_union_matrix) - I ## subtract identity matrix
	replace!(gene_distance_matrix, NaN=>0)
	
	gene_distance_matrix = triu(gene_distance_matrix') + gene_distance_matrix
	genome_distance_matrix = triu(genome_distance_matrix') + genome_distance_matrix
	
	row_hclust = hclust(genome_distance_matrix, branchorder = :optimal)
	col_hclust = hclust(gene_distance_matrix, branchorder = :optimal)
	
	sorted_accessions_seq_matrix = accessions_seq_matrix[row_hclust.order, col_hclust.order]
	
	return sorted_accessions_seq_matrix
	
end

# ╔═╡ ed6bff37-2f4f-4e17-9209-24ba7c900165
accessions_seq_matrix = ProteinToAccessionMatrix(HGT_dup_protseq_to_accession)

# ╔═╡ 306ad03e-6ef6-4964-8ab1-9c51a02c2319
md""" 

after using MinHashing, calculate 1 - Jaccard similarity to get a distance matrix, and 
use this to cluster the rows and then cluster the columns using hierarchical clustering.

Then, use the order of rows and the order of columns returned by hierarchical clustering to sort the rows and columns. Then, visualize the original matrix.


"""

# ╔═╡ c0aaa528-b72c-40fd-a76d-d1aef3a41fa9
begin
	
	HGT_dup_accessions_seq_matrix = SortAccessionSeqMatrix(HGT_dup_protseq_to_accession, HGT_dup_accession_to_protseq) |>
	sparse
	
	spy(HGT_dup_accessions_seq_matrix)
	
end

# ╔═╡ 23217fb2-5a00-4cdc-8191-8b0d4a23e0e0
begin
	
	genomes_sharing_dups_matrix = HGT_dup_accessions_seq_matrix * 	HGT_dup_accessions_seq_matrix'
	
	spy(genomes_sharing_dups_matrix)
	
end

# ╔═╡ 454af6f7-0751-4828-a659-224da4793ec0
## This regex is used in multiple blocks of code.
IS_regex = r"IS|transposon|Transposase|transposase|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid"

# ╔═╡ 531d4292-596d-47c5-a1cb-fb89e01eefca
begin
	## make a dictionary of protein sequence to Annotation_Accession.
	## NOTE: This works on HGT_dup_protein_df by design.
	only_MGE_protseq_to_accession = Dict()
	
	## filter for sequences that are present in multiple genomes,
	## by considering HGT_dup_protein_df.
	for row in eachrow(HGT_dup_protein_df)
		accession = row.Annotation_Accession
		sequence = row.sequence
		annotation = row.product
		if occursin(IS_regex, annotation) ## only consider if matches IS regex.
			if haskey(only_MGE_protseq_to_accession, sequence)
				push!(only_MGE_protseq_to_accession[sequence], accession)
			else
				only_MGE_protseq_to_accession[sequence] = [accession]
			end
		end
	end
	
end

# ╔═╡ d37d69d2-9d9c-40aa-bf16-c7e4eecda90b
only_MGE_HGT_dup_protein_df = filter(row -> occursin(IS_regex, row.product),  HGT_dup_protein_df)

# ╔═╡ 630ee5a6-e1be-4b74-acb3-b13640119e22
## make a dictionary of Annotation_Accession to protein sequence.
only_MGE_HGT_dup_accession_to_protseq = AccessionToProtein(only_MGE_HGT_dup_protein_df)

# ╔═╡ a3aa8cf4-bbb0-42c6-af0f-01ce063aa487
begin	
	
	sparse_only_MGE_accessions_seq_matrix = SortAccessionSeqMatrix(only_MGE_protseq_to_accession, only_MGE_HGT_dup_accession_to_protseq) |>
	sparse
	
	spy(sparse_only_MGE_accessions_seq_matrix)
	
end

# ╔═╡ 39f04c76-c9db-4c76-80db-43ac22e02909
begin
	
	only_MGE_genomes_sharing_dups_matrix = sparse_only_MGE_accessions_seq_matrix * sparse_only_MGE_accessions_seq_matrix'
	
	spy(only_MGE_genomes_sharing_dups_matrix)
	
end

# ╔═╡ 33b1cd3c-b926-457d-b2ab-0506c128446a
begin
	
	only_MGE_shared_gene_dups_matrix = sparse_only_MGE_accessions_seq_matrix' * sparse_only_MGE_accessions_seq_matrix
	
	spy(only_MGE_shared_gene_dups_matrix)
	
end

# ╔═╡ f04c8c61-c9e7-4671-93b8-649ba0b39bb0
md""" let's make a similar matrix plot, but without MGE. """

# ╔═╡ 8651d0e3-e7d0-45e3-b6c5-83bcf179a1d1
begin
	## make a dictionary of protein sequence to Annotation_Accession.
	## NOTE: This works on HGT_dup_protein_df by design.
	no_MGE_protseq_to_accession = Dict()
		
	## filter for sequences that are present in multiple genomes,
	## by considering HGT_dup_protein_df.
	for row in eachrow(HGT_dup_protein_df)
		accession = row.Annotation_Accession
		sequence = row.sequence
		annotation = row.product
		
		if occursin(IS_regex, annotation) ## then skip.
			continue
		end

		if haskey(no_MGE_protseq_to_accession, sequence)
			push!(no_MGE_protseq_to_accession[sequence], accession)
		else
			no_MGE_protseq_to_accession[sequence] = [accession]
		end
	end
	
end

# ╔═╡ 0ed7ead1-b66f-4a75-bdce-97ff3fbe74a7
no_MGE_HGT_dup_protein_df = filter(row -> !occursin(IS_regex, row.product),  HGT_dup_protein_df)

# ╔═╡ fabb3528-65ee-4ee7-95e1-aa1528043553
## make a dictionary of Annotation_Accession to protein sequence.
no_MGE_HGT_dup_accession_to_protseq = AccessionToProtein(no_MGE_HGT_dup_protein_df)

# ╔═╡ f2855669-0695-45ed-9d8c-633c415dc458
begin

	sparse_no_MGE_accessions_seq_matrix = SortAccessionSeqMatrix(no_MGE_protseq_to_accession,
		no_MGE_HGT_dup_accession_to_protseq) |>
	sparse
	
	spy(sparse_no_MGE_accessions_seq_matrix)
	
end

# ╔═╡ b40bd7bc-597b-4f49-bdc7-eda138cb1d12
begin
	
	no_MGE_genomes_sharing_dups_matrix = sparse_no_MGE_accessions_seq_matrix * sparse_no_MGE_accessions_seq_matrix'
	
	spy(no_MGE_genomes_sharing_dups_matrix)
	
end

# ╔═╡ f5accd22-4b48-46ea-9be3-50b97719618f
begin
	
	no_MGE_shared_gene_dups_matrix = sparse_no_MGE_accessions_seq_matrix' * sparse_no_MGE_accessions_seq_matrix
	
	spy(no_MGE_shared_gene_dups_matrix)
	
end

# ╔═╡ d9a9ccc0-c5bf-4a07-8e22-f5bfe7ea149d
md""" let's make a similar matrix plot, but only for ARGs. """

# ╔═╡ 4de58977-1b4d-4d50-93e9-df92c87075b8
## This regex is used in multiple blocks of code.
ARG_regex = r"lactamase|chloramphenicol|quinolone|antibiotic resistance|tetracycline|VanZ"

# ╔═╡ 3ad92137-7f06-4305-99b7-07b4fcdc9fe8
begin
	## make a dictionary of protein sequence to Annotation_Accession.
	## NOTE: This works on HGT_dup_protein_df by design.
	ARG_protseq_to_accession = Dict()
	
	## filter for sequences that are present in multiple genomes,
	## by considering HGT_dup_protein_df.
	for row in eachrow(HGT_dup_protein_df)
		accession = row.Annotation_Accession
		sequence = row.sequence
		annotation = row.product

		
		if !occursin(ARG_regex, annotation) ## then skip.
			continue
		end

		if haskey(ARG_protseq_to_accession, sequence)
			push!(ARG_protseq_to_accession[sequence], accession)
		else
			ARG_protseq_to_accession[sequence] = [accession]
		end
	end
	
end

# ╔═╡ 651429ba-342f-41d8-be01-2f027e3105eb
ARG_HGT_dup_protein_df = filter(row -> occursin(ARG_regex, row.product),  HGT_dup_protein_df)

# ╔═╡ 785a9f69-bfe9-4941-af60-fe4bd0cac639
## make a dictionary of Annotation_Accession to protein sequence.
ARG_HGT_dup_accession_to_protseq = AccessionToProtein(ARG_HGT_dup_protein_df)

# ╔═╡ 4d5f0b2a-8012-468c-9817-d4301ae13841
begin
	
	sparse_ARG_accessions_seq_matrix = SortAccessionSeqMatrix(ARG_protseq_to_accession, ARG_HGT_dup_accession_to_protseq) |> sparse
	
	spy(sparse_ARG_accessions_seq_matrix)
	
end

# ╔═╡ aaf162c4-6e58-4ed4-9361-0b40d8a9a3d2
begin
	
	genomes_sharing_dup_ARGs_matrix = sparse_ARG_accessions_seq_matrix * sparse_ARG_accessions_seq_matrix'
	
	spy(genomes_sharing_dup_ARGs_matrix)
	
end

# ╔═╡ Cell order:
# ╟─73979b34-caea-11eb-3431-8d722c6416b9
# ╟─0c0a7d17-44a2-4ec7-aa83-b1d52920cd2f
# ╠═1ba467c8-5a15-4342-adca-486c970a0c0f
# ╟─cf2d881b-0502-4ee0-8b49-94c33427c263
# ╠═ec356b3d-18a3-4977-8fba-2ca9ee7524a7
# ╟─2afb1c61-50a2-4521-8018-9cc5e9b014b6
# ╠═25f6e9b8-a433-4eb1-92f2-91aad7c65e87
# ╠═d492149c-2fe7-4992-a0bb-3c81f27b1b30
# ╠═463fea4c-ccd1-4b70-9e38-0d8f31a523e1
# ╠═f96f9c72-f05a-4998-8588-555120dc61e9
# ╠═b1509dee-b5a0-4dd2-a1eb-60a82bf21866
# ╟─6b5ef553-0bb4-47d1-b6cf-c6fb898675f3
# ╠═86fc624f-4874-4aa4-acd7-6d70ff9c927f
# ╠═d63763c7-ce0d-4f7d-b7d9-2e9b5a8a5fa9
# ╠═e75548f9-60ba-41c5-ad83-c3a1061b4177
# ╟─8cfe264a-187f-4c63-b5e0-80cdaa284896
# ╠═d495cb9f-02e0-45eb-b5d8-0a091faf9d8e
# ╟─196b9ad6-5251-4534-bb7b-185130589ffb
# ╠═2c068fa2-83a9-466a-be31-3d480e18aa21
# ╠═3736e126-ebb4-4949-83a9-14d7922dc3d3
# ╠═c8c5417f-e7c4-4e68-874f-5ed9d15139da
# ╠═ed6bff37-2f4f-4e17-9209-24ba7c900165
# ╟─306ad03e-6ef6-4964-8ab1-9c51a02c2319
# ╠═c0aaa528-b72c-40fd-a76d-d1aef3a41fa9
# ╠═23217fb2-5a00-4cdc-8191-8b0d4a23e0e0
# ╠═454af6f7-0751-4828-a659-224da4793ec0
# ╠═531d4292-596d-47c5-a1cb-fb89e01eefca
# ╠═d37d69d2-9d9c-40aa-bf16-c7e4eecda90b
# ╠═630ee5a6-e1be-4b74-acb3-b13640119e22
# ╠═a3aa8cf4-bbb0-42c6-af0f-01ce063aa487
# ╠═39f04c76-c9db-4c76-80db-43ac22e02909
# ╠═33b1cd3c-b926-457d-b2ab-0506c128446a
# ╠═f04c8c61-c9e7-4671-93b8-649ba0b39bb0
# ╠═8651d0e3-e7d0-45e3-b6c5-83bcf179a1d1
# ╠═0ed7ead1-b66f-4a75-bdce-97ff3fbe74a7
# ╠═fabb3528-65ee-4ee7-95e1-aa1528043553
# ╠═f2855669-0695-45ed-9d8c-633c415dc458
# ╠═b40bd7bc-597b-4f49-bdc7-eda138cb1d12
# ╠═f5accd22-4b48-46ea-9be3-50b97719618f
# ╟─d9a9ccc0-c5bf-4a07-8e22-f5bfe7ea149d
# ╠═4de58977-1b4d-4d50-93e9-df92c87075b8
# ╠═3ad92137-7f06-4305-99b7-07b4fcdc9fe8
# ╠═651429ba-342f-41d8-be01-2f027e3105eb
# ╠═785a9f69-bfe9-4941-af60-fe4bd0cac639
# ╠═4d5f0b2a-8012-468c-9817-d4301ae13841
# ╠═aaf162c4-6e58-4ed4-9361-0b40d8a9a3d2
