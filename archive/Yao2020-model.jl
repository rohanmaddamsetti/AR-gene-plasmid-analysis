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
	using DataFrames
	using CSV
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

# **Yao et al. 2020**: _Model of Transposon dynamics_

Rohan Maddamsetti and Yi Yao.

Julia implementation by Rohan Maddamsetti, ODE model by Yi Yao.
"""

# ╔═╡ 0497b050-1232-11eb-176c-698843f3368b
md""" 

## **Model description**
We constructed a kinetic model consisting of two non-dimensionalized ordinary differential equations (ODEs) to describe the effect of antibiotics on the population dynamics of transposons moving back and forth from chromosomes to plasmids.


$\frac{dS_c}{dt} = (\mu_c - \eta)S_c$
$\frac{dS_p}{dt} = \mu_p S_p + \eta S_c$

Here, $S_c$ represents the population that carries the transposons on the chromosome only, while $S_p$ population carries the transposons on the plasmids (no matter whether it carries transposon on the chromosome or not). $S_c$ and $S_p$ grow in a logistic manner, and we only simulate the exponential growth phase of the two populations. $\mu_c$ and $\mu_p$ are the growth rates of the populations with chromosomal and plasmid-encoded transposons, respectively. $\eta$ is the translocation rate of the transposon from the chromosome to the plasmids. We assume that the probability that more than one plasmid-based transposon (i.e. two or more simultaneous loss events) is lost in any unit of time is negligible compared to the rate of chromosome-to-plasmid transposition. We ignore $S_p$ to $S_c$ transpositions, because for multi-copy plasmids, a single transposition to the chromosome does not affect the remaining plasmids-- we assume the $S_p$ population cannot be transposed back to the $S_c$ population.

To define the criterion when $S_c$ or $S_p$ population will be dominant in the community, we need to compare $\frac{S_p(t)}{S_c(t)}$ with $\frac{S_p(0)}{S_c(0)}$, where $S_p(t)$ and $S_c(t)$ are the population sizes of each population when the time is equal to $t$, and $S_p(0)$ and $S_c(0)$ are the population sizes of each population when the time is equal to 0.

The previous ODEs can be solved as 

$S_c(t) = S_c(0) e^{(\mu_c - \eta)t}$
$S_p(t) = (S_p(0) + \frac{\eta S_c(0)}{\mu_p - \mu_c + \eta}) e^{\mu_pt} - \frac{\eta S_c(0)}{\mu_p - \mu_c + \eta} e^{(\mu_c - \eta)t}$

If we assume the growth rates of both populations are much greater than the transposition rate, such that $\mu_c >> \eta$ and $\mu_p >> \eta$ then the criterion of the population distribution is: $S_p$ will be enriched if $\mu_c < \mu_p$.

We use Hill equations to model the effect of antibiotics on growth rates. $\mu_c(0)$ and $\mu_p(0)$ are the maximum growth rates of the chromosomal-encoded transposon population and the plasmid-encoded transposon population, respectively. Both strains have their half-maximal effective concentrations of antibiotics $K_c$ or $K_p$, and Hill coefficient $n_p$ and $n_c$. For simplicity, we assume that $n_p = n_c$, so we just use $n$ for both. $[A]$ represents the antibiotic concentration. 

The growth rates under different antibiotic concentrations can be described as:

$\mu_c([A]) = \mu_c(0)\frac{K_c^{n}}{K_c^{n} + [A]^{n}}$
$\mu_p([A]) = \mu_p(0)\frac{K_p^{n}}{K_p^{n} + [A]^{n}}$

We assume that the transposon is costly in the host cells in the absence of antibiotics; that is, $\mu_c(0) > \mu_p(0)$. Furthermore, we assume that the strain with transposons on plasmids will produce a higher yield of resistance proteins per cell, thus causing higher antibiotics tolerance, that is, $K_c < K_p$. If the environment is free of antibiotics, then $\mu_c(0) > \mu_p(0)$, and $S_c$ will always be the dominant population. To guarantee that the $S_p$ population will be dominant when the antibiotic concentration is high enough, we choose the parameters in the simulation to satisfy that $\mu_p$ will be greater than $\mu_c$ when the $[A]$ is high enough.

To predict the fractions of $S_c$ and $S_p$ in different communities with different antibiotic concentrations, we use the same exponential growth model. The community starts with $S_c$ only in the beginning, and the fractions of $S_p$ in the community at $t = 100$ are calculated as the plasmid-based transposons fraction.
"""

# ╔═╡ b30a9116-13c2-11eb-1e53-1dc3847f4384
md""" Our code follows this tutorial: 
[https://diffeq.sciml.ai/stable/tutorials/ode_example](https://diffeq.sciml.ai/stable/tutorials/ode_example)
"""

# ╔═╡ 22eba0f2-13c2-11eb-00d3-898eb79bd25d
function calc_μc(μc₀, Kc, n, A)
	μc₀ * (Kc^n)/(Kc^n + A^n)
end

# ╔═╡ 2dd79174-13c2-11eb-1716-6b1949b4632c
function calc_μp(μp₀, Kp, n, A)
		μp₀ * (Kp^n)/(Kp^n + A^n)
end

# ╔═╡ 2d680c50-13c2-11eb-151b-99bacb19999c
function transposon_dynamics!(du,u,p,t)	
		Sc, Sp = u
		μc₀, μp₀, Kc, Kp, n, A, η = p
		
		μc = calc_μc(μc₀, Kc, n, A)
		μp = calc_μp(μp₀, Kp, n, A)
		
    	du[1] = dSc = (μc - η) * Sc
    	du[2] = dSp = μp * Sp + η * Sc
	end

# ╔═╡ 27d302fe-1237-11eb-0166-1bf9048405e7
begin	
	## initial conditions
	Sc₀ = 0.1
	Sp₀ = 1.0000000000000001E-28
	u₀ = [Sc₀, Sp₀]
	tspan = (0.0,100.0) ## time interval
end

# ╔═╡ 1d2eda16-1963-11eb-21a3-8fdaac79fc97
md""" AntibioticConcentration slider"""

# ╔═╡ 29a14964-1963-11eb-138a-df727f736284
@bind AntibioticConcentration Slider(0:0.1:10,show_value=true)

# ╔═╡ 7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
begin
	## default parameters chosen by Yi.
	μc₀ = 0.7
	μp₀ = 0.5
	Kc = 2.5
	Kp = 10.0 
	n = 2 ## Hill coefficient
	## Antibiotic concentration
	A = AntibioticConcentration ## default value = 5
	η = 1.0E-5 ## rate of transposition from chromosome to plasmid
	
	my_parameters = [μc₀, μp₀, Kc, Kp, n, A, η]
end

# ╔═╡ 80e2e102-13cc-11eb-1619-f999c9e0612a
calc_μc(μc₀, Kc, n, A)

# ╔═╡ 4e32b83a-13c8-11eb-390d-0d3502d6f88e
calc_μp(μp₀, Kp, n, A)

# ╔═╡ cecbfcae-1238-11eb-0353-3905b2919507
prob = ODEProblem(transposon_dynamics!, u₀, tspan, my_parameters)

# ╔═╡ 69daf25e-124b-11eb-1fd1-7bb52f61b420
sol = solve(prob)

# ╔═╡ fa177622-124c-11eb-28e1-d99fe7c076a0
plot(sol,linewidth=2,xaxis="t",label=["test"],layout=(2,1))

# ╔═╡ a384ca20-124d-11eb-0bea-87d588481f08
begin	
	tvec = sol.t
	my_Sc = [x[1] for x in sol.u]
	my_Sp = [x[2] for x in sol.u]
	fraction_on_plasmid_vec = my_Sp./(my_Sp + my_Sc)
end

# ╔═╡ 9695e99c-124e-11eb-06bb-0585c7602429
begin
	my_plot = plot(tvec, fraction_on_plasmid_vec,linewidth=2,xaxis="t",label="fraction on plasmid")
	my_plot
end

# ╔═╡ ba7b560c-13c1-11eb-108f-93adf5ae12b3
md""" for robustness, let's re-run these simulations with random starting parameters. Each set of parameters obeys the constraint that 

$\mu_p > \mu_c$ 

at the given antibiotic concentration $[A]$, which should guarantee that the $S_p$ population will be dominant.
"""

# ╔═╡ a278af92-13c5-11eb-06d1-4df33152bd36
function sample_parameter_vec()
	μc₀ = rand() ## between 0 and 1
	μp₀ = rand() ## between 0 and 1
	Kc = 20 * rand() ## between 0 and 20
	Kp = 20 * rand() ## between 0 and 20
	n = 10 * rand() ## Hill coefficient between 0 and 10
	A = 100 * rand() ## Antibiotic concentration between 0 and 100
	η = 1.0E-5 * rand() ## η is always small

	my_parameters = [μc₀, μp₀, Kc, Kp, n, A, η]
end

# ╔═╡ fdc39818-13c4-11eb-1b4c-278ed2e63eed
function array_of_good_parameter_vecs(N)
	## return an array of N parameter vectors, each satisfying μp > μc.
	vec_array = []
	for i in 1:N
		my_parameters = sample_parameter_vec()
		μc₀, μp₀, Kc, Kp, n, A, η = my_parameters
		μc = calc_μc(μc₀, Kc, n, A)
		μp = calc_μp(μp₀, Kp, n, A)
		
		while μc > μp ## resample if chromosome grows faster than plasmid at this [A]
			my_parameters = sample_parameter_vec()
			μc₀, μp₀, Kc, Kp, n, A, η = my_parameters
			μc = calc_μc(μc₀, Kc, n, A)
			μp = calc_μp(μp₀, Kp, n, A)
		end
		append!(vec_array,[my_parameters])
	end
	return vec_array
end

# ╔═╡ 418230ec-1b17-11eb-0076-2536d542f89b
parameter_vectors = array_of_good_parameter_vecs(500)

# ╔═╡ a79aef24-13ca-11eb-0b07-7face8ad3904
begin
	plot2 = plot()
	for param_vec in parameter_vectors
		prob2 = ODEProblem(transposon_dynamics!, u₀, tspan, param_vec)
		sol2 = solve(prob2)
		tvec = sol2.t
		my_Sc = [x[1] for x in sol2.u]
		my_Sp = [x[2] for x in sol2.u]
		fraction_on_plasmid_vec = my_Sp./(my_Sp + my_Sc)
		plot!(plot2, tvec, fraction_on_plasmid_vec,legend = false)
	end
	xlabel!(plot2, "Time")
	ylabel!(plot2, "Frequency of plasmid-borne ARG")
	savefig(plot2,"../results/Yao2020/random-trajectories.pdf")
	plot2
end

# ╔═╡ 90705bbe-3022-11eb-0ed5-19dbdd500308
let
	plot3 = plot(yaxis=:log)
	for param_vec in parameter_vectors
		prob2 = ODEProblem(transposon_dynamics!, u₀, tspan, param_vec)
		sol2 = solve(prob2)
		tvec = sol2.t
		my_Sc = [x[1] for x in sol2.u]
		my_Sp = [x[2] for x in sol2.u]
		fraction_on_plasmid_vec = my_Sp./(my_Sp + my_Sc)
		plot!(plot3, tvec, fraction_on_plasmid_vec,legend = false)
	end
	xlabel!(plot3, "Time")
	ylabel!(plot3, "log(Frequency of plasmid-borne ARG)")
	savefig(plot3,"../results/Yao2020/random-log-trajectories.pdf")
	plot3
end

# ╔═╡ 2f3d0ebe-3024-11eb-1dbd-b166e07b78fc
begin
	## write random parameters to a CSV
	param_df = DataFrame(mu_c0 = [v[1] for v in parameter_vectors],
					mu_p0 = [v[2] for v in parameter_vectors],
               		Kc = [v[3] for v in parameter_vectors],
					Kp = [v[4] for v in parameter_vectors],
					n = [v[5] for v in parameter_vectors],
					Aconc = [v[6] for v in parameter_vectors],
					eta = [v[7] for v in parameter_vectors]
               		)
	CSV.write("../results/Yao2020/random_parameters.csv", param_df)
end

# ╔═╡ d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
bigbreak = html"<br><br><br><br>";

# ╔═╡ Cell order:
# ╟─49567f8e-09a2-11eb-34c1-bb5c0b642fe8
# ╟─181e156c-0970-11eb-0b77-49b143cc0fc0
# ╠═2dcb18d0-0970-11eb-048a-c1734c6db842
# ╟─0497b050-1232-11eb-176c-698843f3368b
# ╟─b30a9116-13c2-11eb-1e53-1dc3847f4384
# ╠═22eba0f2-13c2-11eb-00d3-898eb79bd25d
# ╠═2dd79174-13c2-11eb-1716-6b1949b4632c
# ╠═2d680c50-13c2-11eb-151b-99bacb19999c
# ╠═7d9153d0-13c2-11eb-1c1e-a7e70aaa9072
# ╠═80e2e102-13cc-11eb-1619-f999c9e0612a
# ╠═4e32b83a-13c8-11eb-390d-0d3502d6f88e
# ╠═27d302fe-1237-11eb-0166-1bf9048405e7
# ╠═cecbfcae-1238-11eb-0353-3905b2919507
# ╠═69daf25e-124b-11eb-1fd1-7bb52f61b420
# ╠═fa177622-124c-11eb-28e1-d99fe7c076a0
# ╟─a384ca20-124d-11eb-0bea-87d588481f08
# ╠═1d2eda16-1963-11eb-21a3-8fdaac79fc97
# ╠═29a14964-1963-11eb-138a-df727f736284
# ╠═9695e99c-124e-11eb-06bb-0585c7602429
# ╟─ba7b560c-13c1-11eb-108f-93adf5ae12b3
# ╠═a278af92-13c5-11eb-06d1-4df33152bd36
# ╠═fdc39818-13c4-11eb-1b4c-278ed2e63eed
# ╠═418230ec-1b17-11eb-0076-2536d542f89b
# ╠═a79aef24-13ca-11eb-0b07-7face8ad3904
# ╠═90705bbe-3022-11eb-0ed5-19dbdd500308
# ╠═2f3d0ebe-3024-11eb-1dbd-b166e07b78fc
# ╟─d5cb6b2c-0a66-11eb-1aff-41d0e502d5e5
