### A Pluto.jl notebook ###
# v0.12.4

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
# **_Toy model of gene duplication dynamics v0.3_**

Rohan Maddamsetti.
"""

# ╔═╡ 9d4909f2-1705-11eb-3bc4-1357471cf06e
md""" 

## **Model description**
"""

# ╔═╡ 6724e69e-16ff-11eb-2b4c-a50d658ddd53
load("../results/AR-gene-duplication/diagrams/toy-model-v0.3.jpg")

# ╔═╡ 0497b050-1232-11eb-176c-698843f3368b
md"""
I built a toy model to illustrate my intuition for why multiple identical copies of a protein sequence in a genome may reveal recent positive selection. The model is shown graphically above, and its equations are below. $x_1$ is the state in the top left corner, and the numbering of states moves roughly clockwise (i.e. $x_2$ is the bottom left state, $x_4$ is the top right state, but $x_5$ is in the bottom right corner).


$\frac{dx_1}{dt} = (f_1 - \phi)x_1 - \mu x_1 - k_{dup} x_1 + k_{rec} x_2 + 0.5 k_{rec} x_3$ 

$\frac{dx_2}{dt} = (f_2 - \phi)x_2 + k_{dup} x_1 - k_{rec} x_2 - 2 \mu x_2$

$\frac{dx_3}{dt} = (f_3 - \phi)x_3 + 2 \mu x_2 - k_{rec} x_3$

$\frac{dx_4}{dt} = (f_4 - \phi)x_4 + \mu x_1 + 0.5 k_{rec} x_3 - k_{dup} x_4 + k_{rec} x_5$

$\frac{dx_5}{dt} = (f_5 - \phi)x_5 + k_{dup} x_4 - k_{rec} x_5$

frequencies add to 1: $x_1 + x_2 + x_3 + x_4 + x_5 = 1$ 

mean fitness: $\phi = \Sigma f_i x_i$

rate of gene duplication: $k_{dup}$

rate of gene loss by recombination: $k_{rec}$

growth rate (fitness) of cell type $i$ : $f_i$. 

I assume that $f_i > 0$. $x_1$ represents the frequency of cells with an antibiotic resistance gene (ARG). $x_2$ represents the frequency of cells with a duplicated ARG. $x_3$ represents the frequency of cell in which one of the two duplicate copies has a beneficial amino acid substitution. $x_4$ represents the frequency of cells in which the mutated cell type has lost the ancestral protein sequence.

Fitness functions:

$f_1 = \frac{K}{K + A}$ where $A$ is antibiotic concentration and $K$ is the concentration of antibiotic that reduces growth by 50%. 

$f_2 = (1-c) \frac{2K}{2K + A}$ where $c$ is the fitness cost of a duplicate gene. 

We assume that $0 < c < 1$, that $A > 0$ and that $K = 1$ for simplicity. We assume that the mutated allele has a benefit $b > 0$. We also assume that the mutated allele is co-dominant (in the sense of Mendelian genetics), such that the heterozygote shows an intermediate phenotype. So,

$f_1 = \frac{1}{1 + A}$

$f_2 = (1-c) \frac{2}{2 + A}$

$f_3 = (1-c) \frac{0.5(1+(1+b))}{0.5(1+(1+b)) + A}$

$f_4 = \frac{(1+b)}{(1+b) + A}$

$f_5 = (1-c) \frac{2(1+b)}{2(1+b) + A}$

I use values from the research literature for $k_{rec}$, $k_{dup}$, and $\mu$.

Tomanek et al. (2020) "Gene amplification as a form of population-level gene expression regulation" gives $k_{rec} = 0.0134(cell^{-1})(generation ^{-1})$ based on their calibration experiments, and $k_{dup} = 0.0001(cell^{-1})(generation ^{-1})$ based on their literature search. Both parameters are listed in their Supplementary Table 2.

I let $\mu = 10^{-7} = 1000(nt) * 10^{-10}(nt^{-1})(cell^{-1})(generation ^{-1})$, where we assume that the ARG is 1000 nucleotides long, and using the point mutation rate reported in Lee H, Popodi E, Tang H, Foster PL. Rate and molecular spectrum of spontaneous mutations in the bacterium Escherichia coli as determined by whole-genome sequencing. Proc Natl Acad Sci U S A. 2012. _Of course, this is a very large upper bound on the rate of beneficial mutations in the antibiotic resistance gene!_

"""

# ╔═╡ b30a9116-13c2-11eb-1e53-1dc3847f4384
md""" My code follows this tutorial: 
[https://diffeq.sciml.ai/stable/tutorials/ode_example](https://diffeq.sciml.ai/stable/tutorials/ode_example)
"""

# ╔═╡ 22eba0f2-13c2-11eb-00d3-898eb79bd25d
function calc_f1(Aₜ)
	1/(1+Aₜ)
end

# ╔═╡ 2dd79174-13c2-11eb-1716-6b1949b4632c
function calc_f2(Aₜ, c)
	(1 - c) * (2/(2 + Aₜ))
end

# ╔═╡ 3761fd48-1934-11eb-357c-4fe040323da7
function calc_f3(Aₜ, c, b)
	Kavg = (1+(1+b))/2
	(1 - c) * (Kavg/(Kavg + Aₜ))
end

# ╔═╡ 38d947a8-1934-11eb-2283-f7f83dc42250
function calc_f4(Aₜ, c, b)
	(1+b)/((1+b) + Aₜ)
end

# ╔═╡ 3985e468-1934-11eb-3517-4f75a7c018b2
function calc_f5(Aₜ, c, b)
	(1 - c) * (2(1+b)/(2(1+b) + Aₜ))
end

# ╔═╡ 42e3448a-170f-11eb-2865-a51c0342b813
function calc_ϕ(u, fvec)
	x1, x2, x3, x4, x5 = u
	f1, f2, f3, f4, f5 = fvec
	ϕ = f1*x1 + f2*x2 + f3*x3 + f4*x4 + f5*x5
	return(ϕ)
end

# ╔═╡ 2d680c50-13c2-11eb-151b-99bacb19999c
function duplication_dynamics!(du,u,p,t)
	x1, x2, x3, x4, x5 = u
	## Note: A is a function of time and not a scalar!
	μ, Krec, Kdup, A, c, b = p

	f1 = calc_f1(A(t))
	f2 = calc_f2(A(t), c)
	f3 = calc_f3(A(t), c, b)
	f4 = calc_f4(A(t), c, b)
	f5 = calc_f5(A(t), c, b)
	f_vec = [f1, f2, f3, f4, f5]

	ϕ = calc_ϕ(u, f_vec)

    du[1] = dx1 = (f1 - ϕ)*x1 - μ*x1 - Kdup*x1 + Krec*x2 + 0.5Krec*x3
    du[2] = dx2 = (f2 - ϕ)*x2 + Kdup*x1 - Krec*x2 - 2μ*x2
	du[3] = dx3 = (f3 - ϕ)*x3 + 2μ*x2 - Krec*x3
	du[4] = dx4 = (f4 - ϕ)*x4 + μ*x1 - Kdup*x4 + Krec*x5 + 0.5Krec*x3 
	du[5] = dx5 = (f5 - ϕ)*x5 + Kdup*x4 - Krec*x5
end

# ╔═╡ 27d302fe-1237-11eb-0166-1bf9048405e7
begin	
	## initial conditions
	x1, x2, x3, x4, x5 = 1,0,0,0,0
	u₀ = [x1, x2, x3, x4,x5]
	## time interval
	tspan = (0.0,1000.0) 
end

# ╔═╡ 567844ce-1936-11eb-002d-3d21a84ec80f
md""" slider for AntibioticConcentration"""

# ╔═╡ 576717f4-1936-11eb-357b-19bd39b2660e
@bind AntibioticConcentration Slider(0:0.1:5,show_value=true)

# ╔═╡ 5830e9ee-1936-11eb-25e1-713d7800be2a
md""" slider for DuplicationCost"""

# ╔═╡ 58e431f0-1936-11eb-3ced-05c917af68c0
@bind DuplicationCost Slider(0:0.1:1,show_value=true)

# ╔═╡ a2f789a0-1935-11eb-3364-6333d30c1979
md""" slider for BeneficialMutationEffect"""

# ╔═╡ adde9eec-1935-11eb-0453-f594e5b0a914
@bind BeneficialMutationEffect Slider(0:0.1:1,show_value=true)

# ╔═╡ 2e9e04c6-193a-11eb-1f5d-cb85d819ce2b
md""" slider for PulseTimeLength"""

# ╔═╡ 3a9ca3c2-193a-11eb-2678-3de07ba44c99
@bind PulseTimeLength Slider(0:10:1000,show_value=true)

# ╔═╡ 7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
begin
	## parameters from the literature
	μ = 10^-7 ## from Lee et al. (2012).
	Krec = 0.0134 ## from Tomanek et al. (2020)
	Kdup = 0.0001 ## from Tomanek et al. (2020)
	## Antibiotic concentration as a function of time.
	A₀ = t->AntibioticConcentration ## constant function.
	Apulse = t->ifelse(t<PulseTimeLength,AntibioticConcentration,0) ## step function.
	## fitness cost of duplication c: can be anywhere between 0 and 1.
	c = DuplicationCost
	## User-defined beneficial mutation effect
	b = BeneficialMutationEffect
	## bundle parameters into a vector.
	antibiotic_treatment = [μ, Krec, Kdup, A₀, c, b]
	pulse_antibiotic_treatment = [μ, Krec, Kdup, Apulse, c, b]
end

# ╔═╡ cecbfcae-1238-11eb-0353-3905b2919507
antibiotic_prob = ODEProblem(duplication_dynamics!, u₀, tspan, antibiotic_treatment)

# ╔═╡ 69daf25e-124b-11eb-1fd1-7bb52f61b420
antibiotic_sol = solve(antibiotic_prob);

# ╔═╡ ec31ed3c-17c7-11eb-0700-6dae1e0fa0ad
pulse_antibiotic_prob = ODEProblem(duplication_dynamics!, u₀, tspan, pulse_antibiotic_treatment)

# ╔═╡ 367b0e5c-17c6-11eb-11ff-25fe3e16ecea
pulse_antibiotic_sol = solve(pulse_antibiotic_prob);

# ╔═╡ fa177622-124c-11eb-28e1-d99fe7c076a0
plot(antibiotic_sol,linewidth=2,xaxis="t")

# ╔═╡ 47619c5e-17c6-11eb-1f77-25bc5f4a5160
plot(pulse_antibiotic_sol,linewidth=2,xaxis="t")

# ╔═╡ d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
bigbreak = html"<br><br><br><br>";

# ╔═╡ Cell order:
# ╟─49567f8e-09a2-11eb-34c1-bb5c0b642fe8
# ╟─181e156c-0970-11eb-0b77-49b143cc0fc0
# ╠═2dcb18d0-0970-11eb-048a-c1734c6db842
# ╟─9d4909f2-1705-11eb-3bc4-1357471cf06e
# ╟─6724e69e-16ff-11eb-2b4c-a50d658ddd53
# ╟─0497b050-1232-11eb-176c-698843f3368b
# ╟─b30a9116-13c2-11eb-1e53-1dc3847f4384
# ╠═22eba0f2-13c2-11eb-00d3-898eb79bd25d
# ╠═2dd79174-13c2-11eb-1716-6b1949b4632c
# ╠═3761fd48-1934-11eb-357c-4fe040323da7
# ╠═38d947a8-1934-11eb-2283-f7f83dc42250
# ╠═3985e468-1934-11eb-3517-4f75a7c018b2
# ╠═42e3448a-170f-11eb-2865-a51c0342b813
# ╠═2d680c50-13c2-11eb-151b-99bacb19999c
# ╠═7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
# ╠═27d302fe-1237-11eb-0166-1bf9048405e7
# ╠═cecbfcae-1238-11eb-0353-3905b2919507
# ╠═ec31ed3c-17c7-11eb-0700-6dae1e0fa0ad
# ╠═69daf25e-124b-11eb-1fd1-7bb52f61b420
# ╠═367b0e5c-17c6-11eb-11ff-25fe3e16ecea
# ╟─567844ce-1936-11eb-002d-3d21a84ec80f
# ╠═576717f4-1936-11eb-357b-19bd39b2660e
# ╟─5830e9ee-1936-11eb-25e1-713d7800be2a
# ╠═58e431f0-1936-11eb-3ced-05c917af68c0
# ╠═a2f789a0-1935-11eb-3364-6333d30c1979
# ╠═adde9eec-1935-11eb-0453-f594e5b0a914
# ╠═2e9e04c6-193a-11eb-1f5d-cb85d819ce2b
# ╠═3a9ca3c2-193a-11eb-2678-3de07ba44c99
# ╠═fa177622-124c-11eb-28e1-d99fe7c076a0
# ╠═47619c5e-17c6-11eb-1f77-25bc5f4a5160
# ╟─d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
