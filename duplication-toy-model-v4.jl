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

model version 4.

Rohan Maddamsetti.
"""

# ╔═╡ 9d4909f2-1705-11eb-3bc4-1357471cf06e
md""" 

## **Model description**
"""

# ╔═╡ 6724e69e-16ff-11eb-2b4c-a50d658ddd53
load("../results/AR-gene-duplication/diagrams/toy-model-v0.4.jpg")

# ╔═╡ 0497b050-1232-11eb-176c-698843f3368b
md"""
I built a toy model to illustrate why multiple identical copies of a protein sequence in a genome may reveal recent positive selection. The model is shown graphically above, and its equations are below.

$\frac{dx_1}{dt} = (f_1 - \phi)x_1 + r x_2 + r x_3$ 

$\frac{dx_2}{dt} = (f_2 - \phi)x_2 - (r + d + \eta)x_2 + r x_4 + r x_5$

$\frac{dx_3}{dt} = (f_3 - \phi)x_3 + \eta x_2 - (r + d)x_3 + r x_5 + r x_6$

$\frac{dx_4}{dt} = (f_4 - \phi)x_4 + d x_2 - (r + \eta)x_4$

$\frac{dx_5}{dt} = (f_5 - \phi)x_5 + \eta x_4 - (2r + \eta)x_5$

$\frac{dx_6}{dt} = (f_6 - \phi)x_6 + d x_3 + \eta x_5 - r x_6$

frequencies add to 1: $\Sigma x_i = 1$ 

mean fitness: $\phi = \Sigma f_i x_i$

rate of gene duplication: $d$

rate of gene loss by recombination: $r$

growth rate (fitness) of cell type $i$ : $f_i$. 

I assume that $f_i > 0$.

Fitness functions:

$f_i = \frac{K}{K + A}$ where $A$ is antibiotic concentration and $K$ is the concentration of antibiotic that reduces growth by 50%. 

$f_j = (1-c) \frac{K}{K + A}$ where $(1-c)$ is the fitness cost of expressing the ARG.

$f_k = (1-c)^2 \frac{K}{K + A}$ where $(1-c)^2$ is the fitness cost of expressing 1 ARG on the plasmid, and an additional copy.

We assume that $0 < c < 1$, that $A > 0$. $K$ varies depending on the configuration of genes on chromosome or plasmid:

$f_1 = \frac{0.1}{0.1 + A}$

$f_2 = \frac{1}{1 + A}$

$f_3 = (1-c) \frac{1.5}{1.5 + A}$

$f_4 = (1-c) \frac{2}{2 + A}$

$f_5 = (1-c)^2 \frac{2.5}{2.5 + A}$

$f_5 = (1-c)^2 \frac{3}{3 + A}$

I use values from the research literature for $r$ and $d$.

Tomanek et al. (2020) "Gene amplification as a form of population-level gene expression regulation" gives $r = 0.0134(cell^{-1})(generation ^{-1})$ based on their calibration experiments, and $d = 0.0001(cell^{-1})(generation ^{-1})$ based on their literature search. Both parameters are listed in their Supplementary Table 2.


"""

# ╔═╡ b30a9116-13c2-11eb-1e53-1dc3847f4384
md""" My code follows this tutorial: 
[https://diffeq.sciml.ai/stable/tutorials/ode_example](https://diffeq.sciml.ai/stable/tutorials/ode_example)
"""

# ╔═╡ 22eba0f2-13c2-11eb-00d3-898eb79bd25d
function calc_f1(Aₜ)
	0.1/(0.1 + Aₜ)
end

# ╔═╡ 2dd79174-13c2-11eb-1716-6b1949b4632c
function calc_f2(Aₜ)
		1/(1 + Aₜ)
end

# ╔═╡ de59b2e6-249a-11eb-08f4-79a001ffd980
function calc_f3(Aₜ, c)
		(1 - c) * 1.5/(1.5 + Aₜ)
end

# ╔═╡ a81f0248-249b-11eb-18bf-5ff78554fe1b
function calc_f4(Aₜ, c)
		(1 - c) * 2/(2 + Aₜ)
end

# ╔═╡ ad33b5a8-249b-11eb-0707-31634e641d07
function calc_f5(Aₜ, c)
		(1 - c)^2 * 2.5/(2.5 + Aₜ)
end

# ╔═╡ ade5699c-249b-11eb-0a96-4b3003350fc6
function calc_f6(Aₜ, c)
		(1 - c)^2 * 3/(3 + Aₜ)
end

# ╔═╡ 42e3448a-170f-11eb-2865-a51c0342b813
function calc_ϕ(u, fitness_vec)
	ϕ = dot(u, fitness_vec)
	return(ϕ)
end

# ╔═╡ 2d680c50-13c2-11eb-151b-99bacb19999c
function dynamics!(du,u,p,t)	
		x1, x2, x3, x4, x5, x6 = u
		## Note: A is a function of time-- and not a scalar!
		r, d, η, A, c = p
		
		f1 = calc_f1(A(t))
		f2 = calc_f2(A(t))
		f3 = calc_f3(A(t), c)
		f4 = calc_f4(A(t), c)
		f5 = calc_f5(A(t), c)
		f6 = calc_f6(A(t), c)
		
		fitness_vec = [f1, f2, f3, f4, f5, f6]
		ϕ = calc_ϕ(u, fitness_vec)
		
    	du[1] = dx1 = (f1 - ϕ)*x1 + r*x2 + r*x3
    	du[2] = dx2 = (f2 - ϕ)*x2 - (r + d + η)*x2 + r*x4 + r*x5
		du[3] = dx3 = (f3 - ϕ)*x3 + η*x2 - (r+d)*x3 + r*x5 + r*x6
		du[4] = dx4 = (f4 - ϕ)*x4 + d*x2 - (r+η)*x4
		du[5] = dx5 = (f5 - ϕ)*x5 + η*x4 - (2r+η)*x5
		du[6] = dx6 = (f6 - ϕ)*x6 + d*x3 + η*x5 - r*x6
end

# ╔═╡ 27d302fe-1237-11eb-0166-1bf9048405e7
begin	
	## initial conditions
	x1, x2, x3, x4, x5, x6 = 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
	u₀ = [x1, x2, x3, x4, x5, x6]
	## time interval
	tspan = (0.0,500.0) 
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

# ╔═╡ 7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
begin 
	## parameters from the literature
	r = 0.0134 ## from Tomanek et al. (2020)
	d = 0.0001 ## from Tomanek et al. (2020)

	## transfer rate.
	η = TransferRate
	
	## Antibiotic concentration as a function of time.
	A₀ = t->AntibioticConcentration ## constant function
	Apulse = t->ifelse(t<200, AntibioticConcentration,0) ## step function
	
	## fitness cost of duplication c: can be anywhere between 0 and 1.
	c = DuplicationCost
	
	## bundle parameters into a vector.
	antibiotic_treatment = [r, d, η, A₀, c]
	pulse_antibiotic_treatment = [r, d, η, Apulse, c]
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

# ╔═╡ bac007fa-1d55-11eb-3075-f5fb5e62ffcf
antibiotic_sol

# ╔═╡ 47619c5e-17c6-11eb-1f77-25bc5f4a5160
begin
p1 = plot([Apulse(t) for t in 1:500], linewidth=1,label="Antibiotic")
p2 = plot(pulse_antibiotic_sol,linewidth=2,xaxis="t")
p3 = plot(p1, p2, layout = (2, 1))
p3
end

# ╔═╡ 3a9451de-1d56-11eb-0011-5df5238d4a71
savefig(p3, "../results/AR-gene-duplication/toy-model-dynamics-v0.4.pdf")

# ╔═╡ 100f8b26-2509-11eb-0166-29830cc4fdf8


# ╔═╡ 0f6b338e-2509-11eb-3fbf-0f09249ff87c


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
# ╟─6724e69e-16ff-11eb-2b4c-a50d658ddd53
# ╟─0497b050-1232-11eb-176c-698843f3368b
# ╟─b30a9116-13c2-11eb-1e53-1dc3847f4384
# ╠═22eba0f2-13c2-11eb-00d3-898eb79bd25d
# ╠═2dd79174-13c2-11eb-1716-6b1949b4632c
# ╠═de59b2e6-249a-11eb-08f4-79a001ffd980
# ╠═a81f0248-249b-11eb-18bf-5ff78554fe1b
# ╠═ad33b5a8-249b-11eb-0707-31634e641d07
# ╠═ade5699c-249b-11eb-0a96-4b3003350fc6
# ╟─42e3448a-170f-11eb-2865-a51c0342b813
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
# ╟─28dfb5a0-249c-11eb-303c-8bd8ebd258e7
# ╠═2e9fded2-249c-11eb-2b78-9be4a0e3e723
# ╠═fa177622-124c-11eb-28e1-d99fe7c076a0
# ╠═bac007fa-1d55-11eb-3075-f5fb5e62ffcf
# ╠═47619c5e-17c6-11eb-1f77-25bc5f4a5160
# ╠═3a9451de-1d56-11eb-0011-5df5238d4a71
# ╠═100f8b26-2509-11eb-0166-29830cc4fdf8
# ╠═0f6b338e-2509-11eb-3fbf-0f09249ff87c
# ╟─c6d3267c-1930-11eb-3164-871a29bfeaad
# ╟─d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
