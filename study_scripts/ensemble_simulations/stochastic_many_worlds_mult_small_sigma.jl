# © 2025 Max Bechthold, John M. Anderies and the IBRI team

using ModelingToolkit, DifferentialEquations, DataFrames, Arrow

"""
This script produces a reduced set of stochastic runs, varying Tg, G0, sigma1 and re1.
It has "large" multiplicative noise (2.5% on average). 
Its results are used to compare to the stochastic runs with small and large noise (figure 6 in the paper).
"""

include("../../src/ibri_std_model.jl")
include("../../src/ibri_analyses.jl")


using .ibri_std_model_mod
using .ibri_analyses_mod



println("Setting up Model...")
ode = ibri_std_model_mod.construct_ode_system()
"""
Equations are:
    D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - r₁₂*P₁*p₁₂ + r₂₁*P₂*p₂₁, #human population dynamics 
    D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + r₁₂*P₁*p₁₂ - r₂₁*P₂*p₂₁, # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - dm₁*K₁, #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - dm₂*K₂, # savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*ct*gmax, # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
"""
ς = 0.025  # small noise level  

noiseeqs = [
    0.0,
    0.0,
    K₁*ς,
    K₂*ς,
    0.0,
    0.0,
    e₁*ς,
    e₂*ς,
]
sde = ibri_std_model_mod.construct_sde_system(ode, noiseeqs)

# lets define the parameters that we want in this study script
tspan=(0.0,610.0); #time span in years. 610 years, such that sims run from 1890 to 2500
u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,
        G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]; #initial conditions
p=[r₁=>0.038, Kp₁=>1.5, #population growth, carrying capacity
r₂=>0.042, Kp₂=>9.7, #population growth, carrying capacity
ml₁₂=>10, ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, #migration preference parameters
r₁₂=>0, r₂₁=>0, #migration rates
α₁=>0.5, α₂=>0.5, #capital factor productivity
a₁=>2.7, a₂=>1.7, #total factor productivity
δ₁=>0.05, δ₂=>0.05, #capital entropic decay rates  
s₁=>0.25, s₂=>0.21, #savings rates (of disposable income)
Cₘ₁=>0.7, Cₘ₂=>0.7, #minimum subsistence consumption
σ₁=>0.03, σ₂=>0.03, dmₓ=> 100, #climate damages on infrastructure
re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>5.0,  #decarb, #carbon intensities
η=>1, #rate at which decarbonization is initiated
u=>0.014, α=>0.1, #earth system carbon uptake and release
G₁=>5, G₀=>5.5, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
];

prob = ibri_std_model_mod.construct_sde_problem(sde, tspan, u0, p)

N_ens = 100
Tg_v = collect(4.2:0.025:6.2)
G0_v = collect(4.2:0.025:6.2)
sigma1_v = [0.025, 0.05, 0.075, 0.1]
re1_v = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15]
N_runs = N_ens * length(Tg_v) * length(G0_v) * length(sigma1_v) * length(re1_v) 

combs = unique(collect(Iterators.product(Tg_v, G0_v, sigma1_v, re1_v)))
further_params = [combs, N_ens, N_runs]
ranges = [Tg_v, G0_v, sigma1_v, re1_v]
ibri_analyses_mod.save_config(".", prob, tspan, u0, p, ranges, further_params)
println("Model setup done!")


println("Start runs...")
sol = ibri_analyses_mod.ensemble_sde_analysis_same_regions_safe_and_just(prob, N_ens, Tg_v, G0_v, sigma1_v, re1_v, [3.5, 4])
println("Runs done!")

println("Saving...")
data_points=[]
for (i, v) in enumerate(combs)
        push!(data_points, (combs[i][1], combs[i][2], combs[i][3], combs[i][4], sol[combs[i]]))
end
df = DataFrame(data_points, [:Tg, :G0, :sigma1, :re1,  :value])
Arrow.write("data/stochastic_many_worlds_mult_small_sigma.arrow", df)
println("Saved!")