### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1ba467c8-5a15-4342-adca-486c970a0c0f
using CSV, DataFrames, Plots, LightGraphs, SparseArrays, Distances, Clustering

# ╔═╡ 73979b34-caea-11eb-3431-8d722c6416b9
md""" ### HGT-network-notebook.jl by Rohan Maddamsetti 

The initial goal of this notebook is to draft matrix and graph visualizations of how identical protein sequences are shared between chromosomes and plasmids, and across strains in my dataset.

I will do analyses that focus just on duplicated genes, as well as singleton genes. I I will then repeat these analyses, focusing on genes associated with mobile genetic elements.

"""

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
begin
	
	## make a matrix as the data structure for the plot.
	## rows are the protein sequences.
	## columns are the Accessions.
	
	accession_set = Set()
	seq_set = Set()
	
	for (seq, accession_list) in HGT_protseq_to_accession
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
	
	for (seq, accession_list) in HGT_protseq_to_accession
		seq_idx = seq_to_index[seq]
		for a in accession_list
			accession_idx = accession_to_index[a]
			accessions_seq_matrix[accession_idx, seq_idx] = 1
		end
	end
	
	
	
end

# ╔═╡ 1455e817-3de0-4906-9b8f-58dc9e3515fa


# ╔═╡ 196b9ad6-5251-4534-bb7b-185130589ffb
md"""  Cluster the rows, and cluster the columns.
 sort and update the seq to index and accession to index mappings
 by getting the permutation of columns and rows returned by the clustering.

 I could potentially speed up these pairwise calculations by parallelizing and/or running on a GPU on the Duke Compute Cluster."""

# ╔═╡ a4a1d61e-b17a-45ad-a878-66000b0cce46
## 620 s to run.
row_dist_matrix = pairwise(Jaccard(), accessions_seq_matrix, dims=1);

# ╔═╡ 653acfac-7436-4e56-855b-9c6924bb2f2e
## 2368 sec to run.
col_dist_matrix = pairwise(Jaccard(), accessions_seq_matrix, dims=2);

# ╔═╡ d598656f-e08b-481a-974c-3d4efab1520c
row_hclust = hclust(row_dist_matrix, linkage=:ward);

# ╔═╡ 51c59c4f-625a-48ec-b14a-fdc833d18b5a
col_hclust = hclust(col_dist_matrix, linkage=:ward);

# ╔═╡ 9faaff9a-67b8-4fa5-9220-0ad1dbf47b3c
sorted_accessions_seq_matrix = accessions_seq_matrix[row_hclust.order,col_hclust.order];

# ╔═╡ ca25b462-b44b-4022-abcd-1c3e53ab4c96
sparse_accessions_seq_matrix = sparse(sorted_accessions_seq_matrix)

# ╔═╡ 8adb4f2e-3d79-4aef-8cde-0074853a2eb7
spy(sparse_accessions_seq_matrix) ## plot the sorted sparse matrix.

# ╔═╡ 23217fb2-5a00-4cdc-8191-8b0d4a23e0e0
genomes_sharing_dups_matrix = sparse_accessions_seq_matrix * sparse_accessions_seq_matrix'

# ╔═╡ 5db94019-2b85-44b9-b8ef-a5ea29eed206
spy(genomes_sharing_dups_matrix)

# ╔═╡ 07084ba8-0dbf-4a3e-abd0-83d4a43f2540
md""" TODO: let's make a similar matrix plot, but only with MGE."""

# ╔═╡ f04c8c61-c9e7-4671-93b8-649ba0b39bb0
md""" let's make a similar matrix plot, but without MGE. """

# ╔═╡ 8651d0e3-e7d0-45e3-b6c5-83bcf179a1d1
begin
	## make a dictionary of protein sequence to Annotation_Accession.
	## NOTE: This works on HGT_dup_protein_df by design.
	no_MGE_protseq_to_accession = Dict()
	
	IS_regex = r"IS|transposon|Transposase|transposase|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid"
	
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

# ╔═╡ 54e24104-f4dd-43a2-b124-9ec1487f06b5
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

# ╔═╡ f2855669-0695-45ed-9d8c-633c415dc458
no_MGE_accessions_seq_matrix = ProteinToAccessionMatrix(no_MGE_protseq_to_accession);

# ╔═╡ 10afcb95-012e-45a9-a966-495505945306
function SortAccessionSeqMatrix(accessions_seq_matrix)
	
	row_dist_matrix = pairwise(Jaccard(), accessions_seq_matrix, dims=1)
	col_dist_matrix = pairwise(Jaccard(), accessions_seq_matrix, dims=2)
	
	row_hclust = hclust(row_dist_matrix, linkage=:ward)
	col_hclust = hclust(col_dist_matrix, linkage=:ward)
	
	sorted_accessions_seq_matrix = accessions_seq_matrix[row_hclust.order,col_hclust.order]
	
	return sorted_accessions_seq_matrix
end

# ╔═╡ 989f05b6-5861-41d1-8b01-5522ea862455
no_MGE_sorted_accessions_seq_matrix = SortAccessionSeqMatrix(no_MGE_accessions_seq_matrix);

# ╔═╡ 8c1086a6-b274-4068-a024-f273f95d2cc3
sparse_no_MGE_accessions_seq_matrix = sparse(no_MGE_sorted_accessions_seq_matrix)

# ╔═╡ c2da23e1-f8c3-4f77-a3d7-47b58e1122f7
spy(sparse_no_MGE_accessions_seq_matrix)

# ╔═╡ b40bd7bc-597b-4f49-bdc7-eda138cb1d12
no_MGE_genomes_sharing_dups_matrix = sparse_no_MGE_accessions_seq_matrix * sparse_no_MGE_accessions_seq_matrix'

# ╔═╡ f1f6396b-6f65-4a57-9983-6821d5dc3b53
spy(no_MGE_genomes_sharing_dups_matrix)

# ╔═╡ f5accd22-4b48-46ea-9be3-50b97719618f
no_MGE_shared_gene_dups_matrix = sparse_no_MGE_accessions_seq_matrix' * sparse_no_MGE_accessions_seq_matrix

# ╔═╡ c301be78-61cf-41ef-91b9-c4d9f5133e3a
spy(no_MGE_shared_gene_dups_matrix)

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
ARG_accessions_seq_matrix = ProteinToAccessionMatrix(ARG_protseq_to_accession);

# ╔═╡ ef7f5f70-f294-43b1-8956-e4eadcdc8365
ARG_sorted_accessions_seq_matrix = SortAccessionSeqMatrix(ARG_accessions_seq_matrix);

# ╔═╡ a1c49718-628d-4fea-a66a-86b72bc16ff6
sparse_ARG_accessions_seq_matrix = sparse(ARG_sorted_accessions_seq_matrix)

# ╔═╡ 181c248a-98e2-40d7-95f7-c1b2017174f8
spy(sparse_ARG_accessions_seq_matrix)

# ╔═╡ aaf162c4-6e58-4ed4-9361-0b40d8a9a3d2
genomes_sharing_dup_ARGs_matrix = sparse_ARG_accessions_seq_matrix * sparse_ARG_accessions_seq_matrix'

# ╔═╡ 5efec77a-235d-4dcc-98ca-25d53ccb67e1
spy(genomes_sharing_dup_ARGs_matrix)

# ╔═╡ 1a3a4c97-a28c-4749-b500-93a2b521e25b
begin
	
	## make the appropriate data structure for the network plot.	
	## plot the network.
	
end

# ╔═╡ Cell order:
# ╟─73979b34-caea-11eb-3431-8d722c6416b9
# ╠═1ba467c8-5a15-4342-adca-486c970a0c0f
# ╟─cf2d881b-0502-4ee0-8b49-94c33427c263
# ╠═ec356b3d-18a3-4977-8fba-2ca9ee7524a7
# ╠═2afb1c61-50a2-4521-8018-9cc5e9b014b6
# ╠═b1509dee-b5a0-4dd2-a1eb-60a82bf21866
# ╠═121b2a92-3867-4ec8-a507-01b6451ed33e
# ╟─6b5ef553-0bb4-47d1-b6cf-c6fb898675f3
# ╠═86fc624f-4874-4aa4-acd7-6d70ff9c927f
# ╠═e75548f9-60ba-41c5-ad83-c3a1061b4177
# ╟─8cfe264a-187f-4c63-b5e0-80cdaa284896
# ╠═d495cb9f-02e0-45eb-b5d8-0a091faf9d8e
# ╠═1455e817-3de0-4906-9b8f-58dc9e3515fa
# ╟─196b9ad6-5251-4534-bb7b-185130589ffb
# ╠═a4a1d61e-b17a-45ad-a878-66000b0cce46
# ╠═653acfac-7436-4e56-855b-9c6924bb2f2e
# ╠═d598656f-e08b-481a-974c-3d4efab1520c
# ╠═51c59c4f-625a-48ec-b14a-fdc833d18b5a
# ╠═9faaff9a-67b8-4fa5-9220-0ad1dbf47b3c
# ╠═ca25b462-b44b-4022-abcd-1c3e53ab4c96
# ╠═8adb4f2e-3d79-4aef-8cde-0074853a2eb7
# ╠═23217fb2-5a00-4cdc-8191-8b0d4a23e0e0
# ╠═5db94019-2b85-44b9-b8ef-a5ea29eed206
# ╠═07084ba8-0dbf-4a3e-abd0-83d4a43f2540
# ╠═f04c8c61-c9e7-4671-93b8-649ba0b39bb0
# ╠═8651d0e3-e7d0-45e3-b6c5-83bcf179a1d1
# ╠═0ed7ead1-b66f-4a75-bdce-97ff3fbe74a7
# ╠═54e24104-f4dd-43a2-b124-9ec1487f06b5
# ╠═f2855669-0695-45ed-9d8c-633c415dc458
# ╠═10afcb95-012e-45a9-a966-495505945306
# ╠═989f05b6-5861-41d1-8b01-5522ea862455
# ╠═8c1086a6-b274-4068-a024-f273f95d2cc3
# ╠═c2da23e1-f8c3-4f77-a3d7-47b58e1122f7
# ╠═b40bd7bc-597b-4f49-bdc7-eda138cb1d12
# ╠═f1f6396b-6f65-4a57-9983-6821d5dc3b53
# ╠═f5accd22-4b48-46ea-9be3-50b97719618f
# ╠═c301be78-61cf-41ef-91b9-c4d9f5133e3a
# ╠═d9a9ccc0-c5bf-4a07-8e22-f5bfe7ea149d
# ╠═3ad92137-7f06-4305-99b7-07b4fcdc9fe8
# ╠═4d5f0b2a-8012-468c-9817-d4301ae13841
# ╠═ef7f5f70-f294-43b1-8956-e4eadcdc8365
# ╠═a1c49718-628d-4fea-a66a-86b72bc16ff6
# ╠═181c248a-98e2-40d7-95f7-c1b2017174f8
# ╠═aaf162c4-6e58-4ed4-9361-0b40d8a9a3d2
# ╠═5efec77a-235d-4dcc-98ca-25d53ccb67e1
# ╠═1a3a4c97-a28c-4749-b500-93a2b521e25b
