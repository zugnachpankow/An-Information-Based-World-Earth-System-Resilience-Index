using ModelingToolkit, DifferentialEquations, Plots, JLD2, DataFrames, Arrow, Random

"""Example code that plots some trajectories of the IBRI standard model."""

# first we include and load the src code needed
include("../src/ibri_std_model.jl")
include("../src/ibri_analyses.jl")

using .ibri_std_model_mod
using .ibri_analyses_mod

println("Setting up Model...")

# use the helper function to construct the ODE system
ode = ibri_std_model_mod.construct_ode_system()

"""
The ODEs are:
    D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - r₁₂*P₁*p₁₂ + r₂₁*P₂*p₂₁, #human population dynamics 
    D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + r₁₂*P₁*p₁₂ - r₂₁*P₂*p₂₁, # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - dm₁*K₁, #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - dm₂*K₂, # savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*ct*gmax, # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
"""

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
re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.2,  #decarb, #carbon intensities
η=>1, #rate at which decarbonization is initiated
u=>0.014, α=>0.1, #earth system carbon uptake and release
G₁=>5, G₀=>4.6, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
];

# now we set up the ODE problem
prob = ibri_std_model_mod.construct_ode_problem(ode, tspan, u0, p)

# construct arrays for changing params
Tg_v = [3.0, 3.4, 3.8, 4.2, 4.6, 5, 5.4, 5.8, 6.2, 6.6, 7.0] # decarbonization points
G0_v = [6.8] # climate threshold
sigma1_v = [0.03] # climate damages
re1_v = [0.1] # decarbonization rate
combs1 = unique(collect(Iterators.product(Tg_v, G0_v, sigma1_v, re1_v)))
println("Model setup done!")


println("Start runs...")
"""a prob_func is needed to change parameters in each run for the ensemble problem
varies the parameters according to the combinations in combs1

prob: the ODE problem
i: the index of the current run
"""
function prob_func1(prob, i, repeat)
        #rand_Tg = rand(pa)
        #rand_G_0 = rand(Gms)
        #rand_Tg = pa[i]
        #rand_G_0 = Gms[i]
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs1[i][3], σ₂=>combs1[i][3], dmₓ=> 100, re₁=>combs1[i][4], 
        re₂=>combs1[i][4], eb=>0.00004, Tg=>combs1[i][1], 
        η=>1, u=>0.0025, α=>0.1, G₁=>5, G₀=>combs1[i][2], Gₘ=>20]
        # push!(rand_values, [rand_Tg, rand_G_0])
        remake(prob, p = pnew) #change Tg and G_0
end

N_samp = length(combs1)
    
# put solver it into a function, since according to documentation thats faster
function construct_solve(prob, prob_func, N_samp)
        ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
        return solve(ensemble_prob, trajectories = N_samp, saveat=1)
end
    
# solve
sol1 = construct_solve(prob, prob_func1, N_samp)
println("Runs done!")

# plot some example trajectories
pp1=plot() # create plot
for i in 1:N_samp
  plot!(pp1,sol1[i][G*100],sol1[i][Y₁+Y₂], xlabel="Climate Threshold G₀ [ppm]",
  ylabel="World Gross National Income [Trillion \$]", alpha=1.0, legend=false, xlim=[280, 800], ylim=[0, 160])
end
display(pp1) # display

