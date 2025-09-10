using ModelingToolkit, DifferentialEquations, Distributions, Roots, NLsolve, DataFrames, Arrow

"""Here we will distribute the climate threshold more realistically."""

include("../../src/ibri_std_model.jl")
include("../../src/ibri_analyses.jl")


using .ibri_std_model_mod
using .ibri_analyses_mod


println("Setting up Model...")
ode = ibri_std_model_mod.construct_ode_system()

# lets define the parameters that we want in this study script

tspan=(0.0,610.0); #time span in years.

u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,
        G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]; #initial conditions

        p=[r₁=>
0.038, Kp₁=>1.5, #population growth, carrying capacity
r₂=>0.042, Kp₂=>9.7, #population growth, carrying capacity
ml₁₂=>10, ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, #migration preference parameters
r₁₂=>0, r₂₁=>0, #migration rates
α₁=>0.5, α₂=>0.5, #capital factor productivity
a₁=>2.7, a₂=>1.7, #total factor productivity
δ₁=>0.05, δ₂=>0.05, #capital entropic decay rates  
s₁=>0.25, s₂=>0.21, #savings rates (of disposable income)
Cₘ₁=>0.7, Cₘ₂=>0.7, #minimum subsistence consumption
σ₁=>0.03, σ₂=>0.03, dmₓ=> 100, #climate damages on infrastructure
re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.5,  #decarb, #carbon intensities
η=>1, #rate at which decarbonization is initiated
u=>0.014, α=>0.1, #earth system carbon uptake and release
G₁=>5, G₀=>4.5, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
];

prob = ibri_std_model_mod.construct_ode_problem(ode, tspan, u0, p)


# construct arrays for changing params
N_samp = 100
Tg_v = collect(4.0:0.01:6.2)
sigma1_v = collect(0.01:0.01:0.1)
re1_v = collect(0.0:0.005:0.1)

combs = collect(Iterators.product(Tg_v, sigma1_v, re1_v))

function temp_to_ppm(x)
        # ppm given as: 5 = 500ppm
        # assume approx 10ppm = 0.1 degree rise
        # so 1° = 100ppm = 1 in our units -> no translation here
        # baseline 280ppm = 2.8, and for temperature its zero since 
        # thresholds are given in anomaly 
        return x + 2.8
end

# Define the equations based on the conditions
function system!(F, x)
        # F[1] = exp(x[1] + x[2]^2 / 2) - 1.6
        F[1] = exp(x[1] + x[2] * quantile(Normal(), 0.05)) - temp_to_ppm(0.8)
        F[2] = exp(x[1] + x[2] * quantile(Normal(), 0.95)) - temp_to_ppm(3.0)
end

solution = nlsolve(system!, [0.0, 1.0])

# Extract the parameters
μ, σ = solution.zero
log_normal_dist = LogNormal(μ, σ)

rand_G0 = rand(truncated(log_normal_dist, 0.001, 100), N_samp)


further_params = [rand_G0]
ranges=[Tg_v, rand_G0, sigma1_v,re1_v]
# for local machine
# ibri_analyses_mod.save_config("Code/studies/high_resolution_many_worlds", prob, tspan, u0, p, further_params)
# for cluster
ibri_analyses_mod.save_config(".", prob, tspan, u0, p, ranges, further_params)
println("Model setup done!")



println("Start runs...")
sol = ibri_analyses_mod.threshold_stochastic_parameters_analysis_safe_and_just(prob, rand_G0, Tg_v, sigma1_v, re1_v, [3.5, 4])
println("Runs done!")

println("Saving...")
# for local machine
# @save "Code/studies/high_resolution_many_worlds/high_resolution_many_worlds.jld2" sol
# for cluster
data_points=[]
for (i, v) in enumerate(combs)
        push!(data_points, (combs[i][1], combs[i][2], combs[i][3], sol[combs[i]]))
end
df = DataFrame(data_points, [:Tg, :sigma1, :re1, :value])
# @save "compressed_results.jld2" sol
Arrow.write("data/deterministic_many_worlds_threshold_distribution.arrow", df)
println("Saved!")







