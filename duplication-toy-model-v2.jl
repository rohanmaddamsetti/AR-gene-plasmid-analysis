### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 2dcb18d0-0970-11eb-048a-c1734c6db842
begin
	using Plots
	using PlutoUI
	using DifferentialEquations
	using LinearAlgebra
	using Images
end

# ╔═╡ 49567f8e-09a2-11eb-34c1-bb5c0b642fe8
# WARNING FOR OLD PLUTO VERSIONS, DONT DELETE ME

html"""
<script>
const warning = html`
<h2 style="color: #800">Oopsie! You need to update Pluto to the latest version for this homework</h2>
<p>Close Pluto, go to the REPL, and type:
<pre><code>julia> import Pkg
julia> Pkg.update("Pluto")
</code></pre>
`

const super_old = window.version_info == null || window.version_info.pluto == null
if(super_old) {
	return warning
}
const version_str = window.version_info.pluto.substring(1)
const numbers = version_str.split(".").map(Number)
console.log(numbers)

if(numbers[0] > 0 || numbers[1] > 12 || numbers[2] > 1) {
	
} else {
	return warning
}

</script>

"""

# ╔═╡ 181e156c-0970-11eb-0b77-49b143cc0fc0
md"""
# **_Toy model of ARG duplication and transposition dynamics_**

**model version 2**.

Rohan Maddamsetti and Lingchong You.
"""

# ╔═╡ 9d4909f2-1705-11eb-3bc4-1357471cf06e
md""" 

## **Model description**
"""

# ╔═╡ 6724e69e-16ff-11eb-2b4c-a50d658ddd53
load("../results/AR-gene-duplication/diagrams/toy-model-v0.5.2.jpg")

# ╔═╡ 0497b050-1232-11eb-176c-698843f3368b
md"""
I built a toy model to illustrate why multiple identical copies of a protein sequence in a genome may reveal recent positive selection. A diagram of the model is shown above.

There are five subpopulations of bacteria. Each cell contains a chromosome and a multi-copy plasmid. Each chromosome and plasmid may contain an antibiotic resistance gene (ARG). An ARG on a chromosome is shown as a large black bar, and an ARG on the plasmid is shown as a small black bar.

We are interested in the dynamics of the five subpopulations due to growth and mutation (duplication, loss, and transfer dynamics of the ARG). I roughly follow the modeling framework used by Lopatkin et al. (2017) "Persistence and reversal of plasmid-mediated antibiotic resistance", and described in that paper's Supplementary Information, and the model used in "Source–sink plasmid transfer dynamics maintain gene mobility in soil bacterial communities" by Hall et al. (2016) in PNAS.

**Model Assumptions**

*Growth dynamics*. We assume that there is a steady inflow of nutrients and antibiotic, and a steady outflow of depleted media and cells, reflected by a constant dilution rate, $D$. This assumption allows the population to grow continously at a steady-state population size. We normalize the number of cells by the carrying capacity, such that each state variable represents the percentage of carrying capacity that is taken up by the subpopulation-- note that this is *not* the relative frequency of cells in the population, because the total population may be at a steady state that is less than carrying capacity. The growth rate of each subpopulation is modeled by  growth functions $f_i > 0$, that we describe in greater detail below.

*Mutation dynamics.* We define a mutation as a transition from one state to another, due to gene duplication, gene loss through recombination, and transposition. Each event occurs at a constant rate $d$, $r$, and $\eta$, respectively. We assume that transposition occurs through a cut-and-paste mechanism, such that the transposition of an ARG from the chromosome to the plasmid removes the ARG from the plasmid. Since we assume a multi-copy plasmid, a transposition from a plasmid to a chromosome removes the ARG from one plasmid copy, but not the others. This justifies the transition from state $x_3$ to state $x_5$.

These assumptions lead to a system of differential equations of the form:

$\frac{dx_i}{dt} = (f_ix_i + Q_i)(1 - \Sigma x_i) - Dx_i$ where the first term reflects logistic growth at a rate $f_i$ when carrying capacity has not been reached, the second term reflects constant dilution due to a fixed outflow rate, and the third term wraps up all the state transitions (mutation dynamics) due to gene duplication, loss, and transposition. 


**Model Equations**


$\frac{dx_1}{dt} = (f_1x_1 - (d + \eta)x_1 + r x_3 + r x_4)(1 - \Sigma x_i) - Dx_1$


$\frac{dx_2}{dt} = (f_2x_2 + \eta x_1 - (d + \eta)x_2 + r x_4 + r x_5) (1 - \Sigma x_i) - Dx_2$

$\frac{dx_3}{dt} = (f_3x_3 + d x_1 - (r + \eta)x_3) (1 - \Sigma x_i) - Dx_3$

$\frac{dx_4}{dt} = (f_4x_4 + \eta x_2 + \eta x_3 - (2r + \eta)x_4) (1 - \Sigma x_i) - Dx_4$

$\frac{dx_5}{dt} = (f_5x_5 + d x_2 + \eta x_4 - r x_5) (1 - \Sigma x_i) - Dx_5$

**Growth functions**

$f_i = \frac{K_i^n}{K_i^n + A^n}$ where $A$ is antibiotic concentration and $K_i$ is the concentration of antibiotic that reduces growth by 50%, and $n$ is a Hill coefficient.

$f_j = (1-c) \frac{K_j^n}{K_j^n + A^n}$ where $(1-c)$ is the cost of expressing the ARG.

$f_k = (1-c)^2 \frac{K_k^n}{K_k^n + A^n}$ where $(1-c)^2$ is the cost of expressing 1 ARG on the plasmid, and an additional copy elsewhere in the genome.

We assume a Hill cofficient $n = 3$. We also assume that $0 < c < 1$, and that $A > 0$. $K$ varies depending on the configuration of genes on chromosome or plasmid:

$f_1 = \frac{1^3}{1^3 + A^3}$

$f_2 = (1-c) \frac{2.5^3}{2.5^3 + A^3}$

$f_3 = (1-c) \frac{5^3}{5^3 + A^3}$

$f_4 = (1-c)^2 \frac{10^3}{10^3 + A^3}$

$f_5 = (1-c)^2 \frac{20^3}{20^3 + A^3}$

I use values from the research literature for $r$ and $d$.

Tomanek et al. (2020) "Gene amplification as a form of population-level gene expression regulation" gives $r = 0.0134(cell^{-1})(generation ^{-1})$ based on their calibration experiments, and $d = 0.0001(cell^{-1})(generation ^{-1})$ based on their literature search. Both parameters are listed in their Supplementary Table 2.


"""

# ╔═╡ b30a9116-13c2-11eb-1e53-1dc3847f4384
md""" My code follows this tutorial: 
[https://diffeq.sciml.ai/stable/tutorials/ode_example](https://diffeq.sciml.ai/stable/tutorials/ode_example)
"""

# ╔═╡ 2dd79174-13c2-11eb-1716-6b1949b4632c
function calc_fi(Aₜ, k, n)
		k^n/(k^n + Aₜ^n)
end

# ╔═╡ de59b2e6-249a-11eb-08f4-79a001ffd980
function calc_fj(Aₜ, k, n, c)
		(1 - c) * k^n/(k^n + Aₜ^n)
end

# ╔═╡ ad33b5a8-249b-11eb-0707-31634e641d07
function calc_fk(Aₜ, k, n, c)
		(1 - c)^2 * k^n/(k^n + Aₜ^n)
end

# ╔═╡ 2d680c50-13c2-11eb-151b-99bacb19999c
function dynamics!(du, u, p, t)	
		x1, x2, x3, x4, x5 = u
		xtotal = sum(u)

		r, d, η, antibiotic_conc_func, c, D = p
		
		k1, k2, k3, k4, k5 = [1, 2.5, 5.0, 10, 20]
		n = 3 # Hill coefficient
	
		f1 = calc_fi(antibiotic_conc_func(t), k1, n)
		f2 = calc_fj(antibiotic_conc_func(t), k2, n, c)
		f3 = calc_fj(antibiotic_conc_func(t), k3, n, c)
		f4 = calc_fk(antibiotic_conc_func(t), k4, n, c)
		f5 = calc_fk(antibiotic_conc_func(t), k5, n, c)
		
    	du[1] = dx1 = (f1*x1 - (d + η)*x1 + r*x3 + r*x4) * (1 - xtotal) - D*x1 
		du[2] = dx2 = (f2*x2 + η*x1 - (d + η)*x2 + r*x4 + r*x5) * (1 - xtotal) - D*x2
		du[3] = dx3 = (f3*x3 + d*x1 - (r+η)*x3) * (1 - xtotal) - D*x3
		du[4] = dx4 = (f4*x4 + η*x2 + η*x3 - (2r+η)*x4) * (1 - xtotal) - D*x4
		du[5] = dx5 = (f5*x5 + d*x2 + η*x4 - r*x5) * (1 - xtotal) - D*x5
end

# ╔═╡ 27d302fe-1237-11eb-0166-1bf9048405e7
begin	
	## initial conditions
	x1, x2, x3, x4, x5 = 0.1, 0.1, 0.1, 0.1, 0.1
	u₀ = [x1, x2, x3, x4, x5]
	## time interval
	tspan = (0.0,1000.0)
end

# ╔═╡ 6bbb5a30-192f-11eb-0dd8-7300290d17f0
md""" Antibiotic Concentration Slider"""

# ╔═╡ 7229eb0c-192f-11eb-0a00-09121752d3c9
@bind AntibioticConcentration Slider(0:0.01:5, show_value=true)

# ╔═╡ 4e324796-1930-11eb-2899-9fa24abd017a
md""" Duplication Cost Slider"""

# ╔═╡ 4f37b380-1930-11eb-248a-975c9d4c9a2e
@bind DuplicationCost Slider(0:0.01:1, show_value=true)

# ╔═╡ 28dfb5a0-249c-11eb-303c-8bd8ebd258e7
md""" Transfer Rate Slider"""

# ╔═╡ 2e9fded2-249c-11eb-2b78-9be4a0e3e723
@bind TransferRate Slider(0:0.00001:0.0002, show_value=true)

# ╔═╡ 673e010e-283e-11eb-288f-fbe7fc2f99bc
md""" Dilution Rate Slider"""

# ╔═╡ 576e9f6a-283e-11eb-0881-df8d68272a19
@bind DilutionRate Slider(0:0.0001:0.2, show_value=true)

# ╔═╡ 7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
begin 
	## parameters from the literature
	r = 0.0134 ## from Tomanek et al. (2020)
	d = 0.0001 ## from Tomanek et al. (2020)

	## transfer rate.
	η = TransferRate
	
	## Dilution rate.
	D = DilutionRate
	
	## Antibiotic concentration as a function of time.
	##A₀ = t->AntibioticConcentration ## constant function
	A₀ = t->AntibioticConcentration
	Apulse = t->ifelse(t<200, AntibioticConcentration,0) ## step function
	
	## fitness cost of duplication c: can be anywhere between 0 and 1.
	c = DuplicationCost
	
	## This condition disables the sliders, and set values to those in
	## our Dynetica model.
	match_Dynetica = true
	if (match_Dynetica)
		η = 3e-5
		D = 0.05
		c = 0.1
		A₀ = t->4
		Apulse = t->ifelse(t<200, 4,0)
	end
	
	## bundle parameters into a vector.
	antibiotic_treatment = [r, d, η, A₀, c, D]
	pulse_antibiotic_treatment = [r, d, η, Apulse, c, D]
end

# ╔═╡ cecbfcae-1238-11eb-0353-3905b2919507
antibiotic_prob = ODEProblem(dynamics!, u₀, tspan, antibiotic_treatment);	

# ╔═╡ 69daf25e-124b-11eb-1fd1-7bb52f61b420
antibiotic_sol = solve(antibiotic_prob);

# ╔═╡ ec31ed3c-17c7-11eb-0700-6dae1e0fa0ad
pulse_antibiotic_prob = ODEProblem(dynamics!, u₀, tspan, pulse_antibiotic_treatment);

# ╔═╡ 367b0e5c-17c6-11eb-11ff-25fe3e16ecea
pulse_antibiotic_sol = solve(pulse_antibiotic_prob);

# ╔═╡ fa177622-124c-11eb-28e1-d99fe7c076a0
plot(antibiotic_sol,linewidth=2,xaxis="t")

# ╔═╡ 381c3af2-3a52-11eb-1353-57f822f177eb
plot(antibiotic_sol,linewidth=2,xaxis="t",yscale=:log10)

# ╔═╡ bac007fa-1d55-11eb-3075-f5fb5e62ffcf
antibiotic_sol

# ╔═╡ 47619c5e-17c6-11eb-1f77-25bc5f4a5160
begin
p1 = plot([Apulse(t) for t in 1:10000], linewidth=1,label="Antibiotic")
p2 = plot(pulse_antibiotic_sol,linewidth=2,xaxis="t")
p3 = plot(p1, p2, layout = (2, 1))
p3
end

# ╔═╡ 3a9451de-1d56-11eb-0011-5df5238d4a71
savefig(p3, "../results/AR-gene-duplication/toy-model-dynamics-v0.5.pdf")

# ╔═╡ 0eeaa15a-4f71-11eb-193c-a10a3c7a964b
md""" __Duplication Index calculation__

DI  (duplication index) = $(x_3 + x_4 + x_5)/ x_total$, that is, the fraction of the gene that will be duplicated (for a certain duplication cost).

DI ~ duplication cost (c) and antibiotic concentration (which controls relative fitnesses).

I can also consider the following alternative metrics:

$x_2 + x_4 + x_5$ fraction of populations containing the gene in the plasmid

$(x_4+x_5)/(x_3+x_4+x_5)$  fraction of duplicated genes in the plasmid


"""

# ╔═╡ 13def1fa-4f71-11eb-20b3-41fecdd8c073
function DuplicationIndex(sol)
	## fraction of population containing duplicates
	return (sol[end][3] + sol[end][4] + sol[end][5])/sum(sol[end])
end

# ╔═╡ 194c7fa4-4f71-11eb-09cf-4f347b7bd421
function PlasmidIndex(sol)
	## fraction of population containing the gene in the plasmid
	return (sol[end][2] +sol[end][4] + sol[end][5])/sum(sol[end])
end

# ╔═╡ 1d679736-4f71-11eb-113a-bd2569486e09
function DuplicatedPlasmidIndex(sol)
	return (sol[end][4] + sol[end][5])/sum(sol[end])
end

# ╔═╡ 22f819a0-4f71-11eb-03bb-d9e161b7fd27
function ConcAndCostToDuplicationIndex(antibiotic_conc::Float64, my_cost::Float64, fixed_parameters)
	
	r, d, η, D = fixed_parameters
	antibiotic_conc_func = t->antibiotic_conc

	my_parameters = [r, d, η, antibiotic_conc_func, my_cost, D]
	
	my_prob = ODEProblem(dynamics!, u₀, tspan, my_parameters)
	my_sol = solve(my_prob)
	return DuplicationIndex(my_sol)
end

# ╔═╡ 27b28d24-4f71-11eb-1291-0b1c965aa0e0
let
	fixed_parameters = [r, d, η, D]
	
	p = plot()
	for cost in 0:0.1:0.3
	antibiotic_concs = [x for x in 0:0.02:5]
	dup_indices = [ConcAndCostToDuplicationIndex(x, cost, fixed_parameters) for x in antibiotic_concs]
	my_label = "cost = $cost"
	plot!(antibiotic_concs, dup_indices, label=my_label,
			legend=:bottomright,ylabel="Duplication Index",
			xlabel="Antibiotic Concentration")
	end
	savefig(p, "../results/AR-gene-duplication/DI-selection-strength.pdf")
	p
end

# ╔═╡ c6d3267c-1930-11eb-3164-871a29bfeaad
md""" **Acknowledgements**

Thanks to Teng Wang and Yi Yao for helpful comments and feedback."""

# ╔═╡ d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
bigbreak = html"<br><br><br><br>";

# ╔═╡ Cell order:
# ╟─49567f8e-09a2-11eb-34c1-bb5c0b642fe8
# ╟─181e156c-0970-11eb-0b77-49b143cc0fc0
# ╠═2dcb18d0-0970-11eb-048a-c1734c6db842
# ╟─9d4909f2-1705-11eb-3bc4-1357471cf06e
# ╠═6724e69e-16ff-11eb-2b4c-a50d658ddd53
# ╟─0497b050-1232-11eb-176c-698843f3368b
# ╟─b30a9116-13c2-11eb-1e53-1dc3847f4384
# ╠═2dd79174-13c2-11eb-1716-6b1949b4632c
# ╠═de59b2e6-249a-11eb-08f4-79a001ffd980
# ╠═ad33b5a8-249b-11eb-0707-31634e641d07
# ╠═2d680c50-13c2-11eb-151b-99bacb19999c
# ╠═7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
# ╠═27d302fe-1237-11eb-0166-1bf9048405e7
# ╠═cecbfcae-1238-11eb-0353-3905b2919507
# ╠═ec31ed3c-17c7-11eb-0700-6dae1e0fa0ad
# ╠═69daf25e-124b-11eb-1fd1-7bb52f61b420
# ╠═367b0e5c-17c6-11eb-11ff-25fe3e16ecea
# ╟─6bbb5a30-192f-11eb-0dd8-7300290d17f0
# ╟─7229eb0c-192f-11eb-0a00-09121752d3c9
# ╟─4e324796-1930-11eb-2899-9fa24abd017a
# ╟─4f37b380-1930-11eb-248a-975c9d4c9a2e
# ╠═28dfb5a0-249c-11eb-303c-8bd8ebd258e7
# ╠═2e9fded2-249c-11eb-2b78-9be4a0e3e723
# ╟─673e010e-283e-11eb-288f-fbe7fc2f99bc
# ╠═576e9f6a-283e-11eb-0881-df8d68272a19
# ╠═fa177622-124c-11eb-28e1-d99fe7c076a0
# ╠═381c3af2-3a52-11eb-1353-57f822f177eb
# ╠═bac007fa-1d55-11eb-3075-f5fb5e62ffcf
# ╠═47619c5e-17c6-11eb-1f77-25bc5f4a5160
# ╠═3a9451de-1d56-11eb-0011-5df5238d4a71
# ╠═0eeaa15a-4f71-11eb-193c-a10a3c7a964b
# ╠═13def1fa-4f71-11eb-20b3-41fecdd8c073
# ╠═194c7fa4-4f71-11eb-09cf-4f347b7bd421
# ╠═1d679736-4f71-11eb-113a-bd2569486e09
# ╠═22f819a0-4f71-11eb-03bb-d9e161b7fd27
# ╠═27b28d24-4f71-11eb-1291-0b1c965aa0e0
# ╟─c6d3267c-1930-11eb-3164-871a29bfeaad
# ╟─d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
