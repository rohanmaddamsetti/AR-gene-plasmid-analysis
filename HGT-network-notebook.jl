### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 1ba467c8-5a15-4342-adca-486c970a0c0f
using CSV, DataFrames, Plots, LightGraphs, SparseArrays, Distances, Clustering

# ╔═╡ 73979b34-caea-11eb-3431-8d722c6416b9
md""" ### HGT-network-notebook.jl by Rohan Maddamsetti 

The initial goal of this notebook is to draft matrix and graph visualizations of how identical protein sequences are shared across strains in my dataset.

I will do analyses that focus just on duplicated genes, as well as singleton genes. I I will then repeat these analyses, focusing on genes associated with mobile genetic elements.

"""

# ╔═╡ 0c0a7d17-44a2-4ec7-aa83-b1d52920cd2f
md""" IMPORTANT TODO!!! Use MinHash algorithm (Mash software) in order to rapidly estimate distances and cluster genes and strains. I need something that scales!!! """

# ╔═╡ e523bcbb-fc4f-4683-a1df-5d2e93f2cf8a
md""" potential TODOs: """

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

# ╔═╡ b1509dee-b5a0-4dd2-a1eb-60a82bf21866
begin
	## make a dictionary of protein sequence to Annotation_Accession.
	protseq_to_accession = Dict()
	
	for row in eachrow(dup_protein_df)
		accession = row.Annotation_Accession
		sequence = row.sequence

		if haskey(protseq_to_accession, sequence)
			push!(protseq_to_accession[sequence], accession)
		else
			protseq_to_accession[sequence] = [accession]
		end
	end
	
end

# ╔═╡ 121b2a92-3867-4ec8-a507-01b6451ed33e
begin
	## let's filter for sequences that are present in multiple genomes.
	HGT_protseq_to_accession = Dict()
	for (k,v) in protseq_to_accession
		if length(v) > 1
			HGT_protseq_to_accession[k] = v
		end
	end
end

# ╔═╡ 6b5ef553-0bb4-47d1-b6cf-c6fb898675f3
md"""

A good heuristic for relatedness is whether two organisms have identical copies of the tuf gene encoding the Tu elongation factor. If two organisms do not share identical tuf sequences, but they have other identical proteins, then it is likely that those identical proteins stem from HGT rather than vertical descent. The same heuristic can be extended to conserved phylogenetic marker genes, etc.

"""

# ╔═╡ 86fc624f-4874-4aa4-acd7-6d70ff9c927f
HGT_dup_protein_df = filter(row -> haskey(HGT_protseq_to_accession,row.sequence), dup_protein_df);

# ╔═╡ e75548f9-60ba-41c5-ad83-c3a1061b4177
plasmid_HGT_dup_protein_df = filter(row -> row.plasmid_count > 1, HGT_dup_protein_df);

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

 I could potentially speed up these pairwise calculations by parallelizing and/or running on a GPU on the Duke Compute Cluster."""

# ╔═╡ 22d9667e-b208-4840-b1d2-43c6cf755a24
function SortAccessionSeqMatrix(accessions_seq_matrix)
	
	row_dist_matrix = pairwise(CosineDist(), accessions_seq_matrix, dims=1)
	col_dist_matrix = pairwise(CosineDist(), accessions_seq_matrix, dims=2)
	
	row_hclust = hclust(row_dist_matrix, linkage=:ward, branchorder = :optimal)
	col_hclust = hclust(col_dist_matrix, linkage=:ward, branchorder = :optimal)
	
	sorted_accessions_seq_matrix = accessions_seq_matrix[row_hclust.order,col_hclust.order]
	
	return sorted_accessions_seq_matrix
	
end

# ╔═╡ 63bc0a20-26d2-40f3-a2c6-594c597ac1a6
begin
	
	## NOTE: SortAccessionSeqMatrix() crashes the worker process,
	## if run using protseq_to_accession as input.
	sparse_accessions_seq_matrix = HGT_protseq_to_accession |>
	ProteinToAccessionMatrix |>
	SortAccessionSeqMatrix |>
	sparse
	
	spy(sparse_accessions_seq_matrix) ## plot the sorted sparse matrix.
end

# ╔═╡ 23217fb2-5a00-4cdc-8191-8b0d4a23e0e0
begin
	
	genomes_sharing_dups_matrix = sparse_accessions_seq_matrix * 	sparse_accessions_seq_matrix'
	
	spy(genomes_sharing_dups_matrix)
	
end

# ╔═╡ 07084ba8-0dbf-4a3e-abd0-83d4a43f2540
md""" TODO: let's make a similar matrix plot, but only with MGE."""

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

# ╔═╡ a3aa8cf4-bbb0-42c6-af0f-01ce063aa487
begin

	sparse_only_MGE_accessions_seq_matrix = only_MGE_protseq_to_accession |> ProteinToAccessionMatrix |> 
	SortAccessionSeqMatrix |>
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

# ╔═╡ f2855669-0695-45ed-9d8c-633c415dc458
begin

	sparse_no_MGE_accessions_seq_matrix = no_MGE_protseq_to_accession |> 	 ProteinToAccessionMatrix |>
	SortAccessionSeqMatrix |>
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
		
		ARG_regex = r"lactamase|chloramphenicol|quinolone|antibiotic resistance|tetracycline|VanZ"
		
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

# ╔═╡ 4d5f0b2a-8012-468c-9817-d4301ae13841
begin
	
	ARG_accessions_seq_matrix = ProteinToAccessionMatrix(ARG_protseq_to_accession)
	
	ARG_sorted_accessions_seq_matrix = SortAccessionSeqMatrix(ARG_accessions_seq_matrix)
	
	sparse_ARG_accessions_seq_matrix = sparse(ARG_sorted_accessions_seq_matrix)
	spy(sparse_ARG_accessions_seq_matrix)
	
end

# ╔═╡ aaf162c4-6e58-4ed4-9361-0b40d8a9a3d2
begin
	
	genomes_sharing_dup_ARGs_matrix = sparse_ARG_accessions_seq_matrix * sparse_ARG_accessions_seq_matrix'
	
	spy(genomes_sharing_dup_ARGs_matrix)
	
end

# ╔═╡ ea5a9e72-034f-407a-905a-424cc9069293
begin
	
	## for comparison, let's look at singletons.
	
	all_protein_fpath = "../results/AR-gene-duplication/all-proteins.csv"
	## get all proteins, then filter for singletons in-place.
	singleton_df = DataFrame(CSV.File(all_protein_fpath))
	filter!(row -> row.count == 1, singleton_df)

end

# ╔═╡ 1a3a4c97-a28c-4749-b500-93a2b521e25b
begin
	
	## make the appropriate data structure for the network plot.	
	## plot the network.
	
end

# ╔═╡ Cell order:
# ╟─73979b34-caea-11eb-3431-8d722c6416b9
# ╠═0c0a7d17-44a2-4ec7-aa83-b1d52920cd2f
# ╠═e523bcbb-fc4f-4683-a1df-5d2e93f2cf8a
# ╠═1ba467c8-5a15-4342-adca-486c970a0c0f
# ╟─cf2d881b-0502-4ee0-8b49-94c33427c263
# ╠═ec356b3d-18a3-4977-8fba-2ca9ee7524a7
# ╟─2afb1c61-50a2-4521-8018-9cc5e9b014b6
# ╠═b1509dee-b5a0-4dd2-a1eb-60a82bf21866
# ╠═121b2a92-3867-4ec8-a507-01b6451ed33e
# ╟─6b5ef553-0bb4-47d1-b6cf-c6fb898675f3
# ╠═86fc624f-4874-4aa4-acd7-6d70ff9c927f
# ╠═e75548f9-60ba-41c5-ad83-c3a1061b4177
# ╟─8cfe264a-187f-4c63-b5e0-80cdaa284896
# ╠═d495cb9f-02e0-45eb-b5d8-0a091faf9d8e
# ╟─196b9ad6-5251-4534-bb7b-185130589ffb
# ╠═22d9667e-b208-4840-b1d2-43c6cf755a24
# ╠═63bc0a20-26d2-40f3-a2c6-594c597ac1a6
# ╠═23217fb2-5a00-4cdc-8191-8b0d4a23e0e0
# ╠═07084ba8-0dbf-4a3e-abd0-83d4a43f2540
# ╠═454af6f7-0751-4828-a659-224da4793ec0
# ╠═531d4292-596d-47c5-a1cb-fb89e01eefca
# ╠═d37d69d2-9d9c-40aa-bf16-c7e4eecda90b
# ╠═a3aa8cf4-bbb0-42c6-af0f-01ce063aa487
# ╠═39f04c76-c9db-4c76-80db-43ac22e02909
# ╠═33b1cd3c-b926-457d-b2ab-0506c128446a
# ╠═f04c8c61-c9e7-4671-93b8-649ba0b39bb0
# ╠═8651d0e3-e7d0-45e3-b6c5-83bcf179a1d1
# ╠═0ed7ead1-b66f-4a75-bdce-97ff3fbe74a7
# ╠═f2855669-0695-45ed-9d8c-633c415dc458
# ╠═b40bd7bc-597b-4f49-bdc7-eda138cb1d12
# ╠═f5accd22-4b48-46ea-9be3-50b97719618f
# ╟─d9a9ccc0-c5bf-4a07-8e22-f5bfe7ea149d
# ╠═3ad92137-7f06-4305-99b7-07b4fcdc9fe8
# ╠═4d5f0b2a-8012-468c-9817-d4301ae13841
# ╠═aaf162c4-6e58-4ed4-9361-0b40d8a9a3d2
# ╠═ea5a9e72-034f-407a-905a-424cc9069293
# ╠═1a3a4c97-a28c-4749-b500-93a2b521e25b
