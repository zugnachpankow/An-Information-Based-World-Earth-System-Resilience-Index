# © 2025 Max Bechthold, John M. Anderies and the IBRI team

using ModelingToolkit, DifferentialEquations, DataFrames, Arrow

"""
This script produces a reduced set of deterministic runs, varying Tg, G0, sigma1 and re1.
Its results are used to compare to the stochastic runs with small and large noise (figure 6 in the paper).
"""

include("../../src/ibri_std_model.jl")
include("../../src/ibri_analyses.jl")


using .ibri_std_model_mod
using .ibri_analyses_mod


println("Setting up Model...")
ode = ibri_std_model_mod.construct_ode_system()

# lets define the parameters that we want in this study script

tspan=(0.0,610.0); #time span in years. 610 years, such that sims run from 1890 to 2500

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

# set up the problem with the helper function
prob = ibri_std_model_mod.construct_ode_problem(ode, tspan, u0, p)

# construct arrays for changing params
Tg_v = collect(4.2:0.025:6.2)
G0_v = collect(4.2:0.025:6.2)
sigma1_v = [0.025, 0.05, 0.075, 0.1]
re1_v = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15]

# construct combinations of parameters
combs = unique(collect(Iterators.product(Tg_v, G0_v, sigma1_v, re1_v)))
N_samp = length(combs)
println("Model setup done!")

println("Start runs...")
sol = ibri_analyses_mod.many_worlds_analysis_same_regions_safe_and_just(prob, combs, [3.5, 4])
println("Runs done!")

println("Saving...")
data_points=[]
for (i, v) in enumerate(combs)
        push!(data_points, (combs[i][1], combs[i][2], combs[i][3], combs[i][4], sol[combs[i]]))
end
df = DataFrame(data_points, [:Tg, :G0, :sigma1, :re1, :value])
Arrow.write("data/deterministic_many_worlds_compare.arrow", df)
println("Saved!")







