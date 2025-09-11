# © 2025 Max Bechthold, John M. Anderies and the IBRI team

using ModelingToolkit, DifferentialEquations, Plots, JLD2, DataFrames, Arrow, Random, Colors, Plots.PlotMeasures, ColorSchemes, Random

"""This script plots Figure 5 as found in the publication:
It creates a plot showing example trajectories of the stochastic versions of the model,
as well as a combination where the Earth system goes into a drifting state (i.e. looses its capacity for carbon uptake)
upon crossing the climate threshold.
To avoid namespace issues we include and use the models inside a module each.
NOTE: we set a random seed for reproducibility.
"""

# setting a random seed for reproducibility
Random.seed!(1234567)

module StandardModels
using ModelingToolkit, DifferentialEquations, JLD2
include("../../src/ibri_std_model.jl")
using .ibri_std_model_mod

function run_std_model()
    println("Setting up and solving Deterministic Standard Model...")
    # Parameter combinations
    Tg_v1 = [4.6, 5.4]
    G0_v1 = [5.0]
    sigma1_v1 = [0.03]
    re1_v1 = [0.1]
    combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
    N_samp = length(combs1)
    # Construct system and problem
    std_ode = ibri_std_model_mod.construct_ode_system()
    tspan = (0.0, 610.0)
    u0 = [P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]
    p = [r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
        mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7,
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7,
        σ₁=>0.03, σ₂=>0.03, dmₓ=>100, re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.2,
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>4.6, Gₘ=>20]
    std_prob = ibri_std_model_mod.construct_ode_problem(std_ode, tspan, u0, p)
    
    # Function to update parameters per sample
    function std_prob_func(prob, i, repeat)
        pnew = [r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10,
                ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5,
                a₁=>2.7, a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21,
                Cₘ₁=>0.7, Cₘ₂=>0.7, σ₁=>combs1[i][3], σ₂=>combs1[i][3],
                dmₓ=>100, re₁=>combs1[i][4], re₂=>combs1[i][4], eb=>0.00004,
                Tg=>combs1[i][1], η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2],
                Gₘ=>20]
        remake(prob, p=pnew)
    end

    # Solve
    ensemble_prob = EnsembleProblem(std_prob, prob_func=std_prob_func)
    return solve(ensemble_prob, trajectories=N_samp, saveat=1)
end # end of function

function run_stochastic_model()
    println("Setting up and solving Stochastic Standard Model...")
    # Load deterministic standard problem first
    Tg_v1 = [4.6, 5.4]; G0_v1 = [5.0]; sigma1_v1 = [0.03]; re1_v1 = [0.1]
    combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
    N_samp = length(combs1)
    tspan=(0.0,610.0)
    u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]
    p=[r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
    mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7,
    a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7,
    σ₁=>0.03, σ₂=>0.03, dmₓ=>100, re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.2,
    η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>4.6, Gₘ=>20]
    std_ode = ibri_std_model_mod.construct_ode_system()
    stochastic_noiseeqs = [0.0, 0.0, K₁*0.025, K₂*0.025, 0.0, 0.0, e₁*0.025, e₂*0.025]
    sde = ibri_std_model_mod.construct_sde_system(std_ode, stochastic_noiseeqs)
    sde_prob = ibri_std_model_mod.construct_sde_problem(sde, tspan, u0, p)
    
    function std_prob_func(prob, i, repeat)
        pnew = [r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
                mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5,
                a₁=>2.7, a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21,
                Cₘ₁=>0.7, Cₘ₂=>0.7, σ₁=>combs1[i][3], σ₂=>combs1[i][3],
                dmₓ=>100, re₁=>combs1[i][4], re₂=>combs1[i][4], eb=>0.00004,
                Tg=>combs1[i][1], η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2],
                Gₘ=>20]
        remake(prob, p=pnew)
    end

    ensemble_prob = EnsembleProblem(sde_prob, prob_func=std_prob_func)
    return solve(ensemble_prob, trajectories=N_samp, saveat=1)
end # end of function

end # end of module


module DriftModels
using ModelingToolkit, DifferentialEquations, JLD2
include("../../src/ibri_drift_model.jl")
using .ibri_drift_model_mod

function run_std_drift_model()
    println("Setting up and solving Deterministic Drift Model...")
    # Parameter combinations (same as before)
    Tg_v1 = [4.6, 5.4]; G0_v1 = [5.0]; sigma1_v1 = [0.03]; re1_v1 = [0.1]
    combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
    N_samp = length(combs1)
    tspan=(0.0,610.0)
    u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]
    p=[r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
    mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7,
    a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7,
    σ₁=>0.03, σ₂=>0.03, dmₓ=>100, re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.2,
    η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>4.6, Gₘ=>20]
    drift_ode = ibri_drift_model_mod.construct_ode_system()
    drift_prob = ibri_drift_model_mod.construct_ode_problem(drift_ode, tspan, u0, p)
    
    function drift_prob_func(prob, i, repeat)
        pnew = [r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
                mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5,
                a₁=>2.7, a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21,
                Cₘ₁=>0.7, Cₘ₂=>0.7, σ₁=>combs1[i][3], σ₂=>combs1[i][3],
                dmₓ=>100, re₁=>combs1[i][4], re₂=>combs1[i][4], eb=>0.00004,
                Tg=>combs1[i][1], η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2],
                Gₘ=>20]
        remake(prob, p=pnew)
    end
    
    ensemble_prob = EnsembleProblem(drift_prob, prob_func=drift_prob_func)
    return solve(ensemble_prob, trajectories=N_samp, saveat=1)
end # end of function

function run_stochastic_drift_model()
    println("Setting up and solving Stochastic Drift Model...")
    # Parameters and combinations
    Tg_v1 = [4.6, 5.4]; G0_v1 = [5.0]; sigma1_v1 = [0.03]; re1_v1 = [0.1]
    combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
    N_samp = length(combs1)
    tspan=(0.0,610.0)
    u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]
    p=[r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
    mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7,
    a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7,
    σ₁=>0.03, σ₂=>0.03, dmₓ=>100, re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.2,
    η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>4.6, Gₘ=>20]
    drift_ode = ibri_drift_model_mod.construct_ode_system()
    stochastic_noiseeqs = [0.0, 0.0, K₁*0.025, K₂*0.025, 0.0, 0.0, e₁*0.025, e₂*0.025]
    drift_sde = ibri_drift_model_mod.construct_sde_system(drift_ode, stochastic_noiseeqs)
    drift_sde_prob = ibri_drift_model_mod.construct_sde_problem(drift_sde, tspan, u0, p)
    
    function drift_prob_func(prob, i, repeat)
        pnew = [r₁=>0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, ml₂₁=>10,
                mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5,
                a₁=>2.7, a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21,
                Cₘ₁=>0.7, Cₘ₂=>0.7, σ₁=>combs1[i][3], σ₂=>combs1[i][3],
                dmₓ=>100, re₁=>combs1[i][4], re₂=>combs1[i][4], eb=>0.00004,
                Tg=>combs1[i][1], η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2],
                Gₘ=>20]
        remake(prob, p=pnew)
    end
    
    ensemble_prob = EnsembleProblem(drift_sde_prob, prob_func=drift_prob_func)
    return solve(ensemble_prob, trajectories=N_samp, saveat=1)    
end # end of function

end # end of module


module ShockModels
using ModelingToolkit, DifferentialEquations, JLD2
include("../../src/ibri_shock_model.jl")
using .ibri_shock_model_mod

function run_shock_model()
    println("Setting up and solving Stochastic Shock Model...")
    Tg_v1 = [4.6, 5.4]
    G0_v1 = [5.0]
    sigma1_v1 = [0.03]
    re1_v1 = [0.1]
    combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
    N_samp = length(combs1)
    tspan=(0.0,610.0); #time span in years.
    shock_ode = ibri_shock_model_mod.construct_ode_system()
    shock_noiseeqs = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    shock_sde = ibri_shock_model_mod.construct_sde_system(shock_ode, shock_noiseeqs)
    shock_u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,
            G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0, a₁=>2.7, a₂=>1.7]; #initial conditions
    global shock_p=[r₁=>
    0.038, Kp₁=>1.5, #population growth, carrying capacity
    r₂=>0.042, Kp₂=>9.7, #population growth, carrying capacity
    ml₁₂=>10, ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, #migration preference parameters
    r₁₂=>0, r₂₁=>0, #migration rates
    α₁=>0.5, α₂=>0.5, #capital factor productivity at max
    an₁=>2.7, an₂=>1.7, #total factor productivity at max
    γ₁=>0.1, γ₂=>0.1, #relaxation speed towards an₁ and an₂
    δ₁=>0.05, δ₂=>0.05, #capital entropic decay rates  
    s₁=>0.25, s₂=>0.21, #savings rates (of disposable income)
    Cₘ₁=>0.7, Cₘ₂=>0.7, #minimum subsistence consumption
    σ₁=>0.03, σ₂=>0.03, dmₓ=> 100, #climate damages on infrastructure
    re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>20,  #decarb, #carbon intensities
    η=>1, #rate at which decarbonization is initiated
    u=>0.014, α=>0.1, #earth system carbon uptake and release
    G₁=>5, G₀=>20, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
    ];
    shock_prob = ibri_shock_model_mod.construct_sde_problem(shock_sde, tspan, shock_u0, shock_p)

    function shock_prob_func(prob, i, repeat)
            #rand_Tg = rand(pa)
            #rand_G_0 = rand(Gms)
            #rand_Tg = pa[i]
            #rand_G_0 = Gms[i]
            pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
            ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, an₁=>2.7, an₂=>1.7,
            δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
            σ₁=>combs1[i][3], σ₂=>combs1[i][3], dmₓ=> 100, re₁=>combs1[i][4], 
            re₂=>combs1[i][4], eb=>0.00004, Tg=>combs1[i][1], 
            η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2], Gₘ=>20]
            # push!(rand_values, [rand_Tg, rand_G_0])
            remake(prob, p = pnew) #change Tg and G_0
    end

    # build jump problem
    λ₀=0.05
    k=1
    Gₚ=3
    crate1(u, p, t) = ibri_shock_model_mod.expit(u[5], λ₀, k, Gₚ)
    crate2(u, p, t) = ibri_shock_model_mod.expit(u[5], λ₀, k, Gₚ)

    function affect1!(integrator)
        integrator.u[9] -= rand(0:0.01:0.2)*integrator.u[9]
        nothing
    end

    function affect2!(integrator)
        integrator.u[10] -= rand(0:0.01:0.2)*integrator.u[10]
        nothing
    end

    crj1 = ConstantRateJump(crate1, affect1!)
    crj2 = ConstantRateJump(crate2, affect2!)
    shock_jump_prob = JumpProblem(shock_prob, Direct(), crj1, crj2) 
    # put it into a function, since according to documentation thats faster
    # Unified construct_solve using output_func
    function construct_solve(prob, prob_func, N_samp)
        ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
        return solve(ensemble_prob, trajectories = N_samp, saveat=1, SRIW1())
    end

    # Solve the shock jump problem
    return construct_solve(shock_jump_prob, shock_prob_func, N_samp) 
end # end of function

end # end of module


module ShockDriftModels
using ModelingToolkit, DifferentialEquations, JLD2
include("../../src/ibri_shock_drift_model.jl")
using .ibri_shock_drift_model_mod

function run_shock_drift_model()
    println("Setting up and solving Shock Drift Model...")
    Tg_v1 = [4.6, 5.4]
    G0_v1 = [5.0]
    sigma1_v1 = [0.03]
    re1_v1 = [0.1]
    combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
    N_samp = length(combs1)
    tspan=(0.0,610.0); #time span in years.
    shock_drift_ode = ibri_shock_drift_model_mod.construct_ode_system()
    # add noise (empty for shocks)
    shock_drift_noiseeqs = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    shock_drift_sde = ibri_shock_drift_model_mod.construct_sde_system(shock_drift_ode, shock_drift_noiseeqs)
    shock_drift_u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,
            G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0, a₁=>2.7, a₂=>1.7]; #initial conditions
    global shock_drift_p=[r₁=>
    0.038, Kp₁=>1.5, #population growth, carrying capacity
    r₂=>0.042, Kp₂=>9.7, #population growth, carrying capacity
    ml₁₂=>10, ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, #migration preference parameters
    r₁₂=>0, r₂₁=>0, #migration rates
    α₁=>0.5, α₂=>0.5, #capital factor productivity at max
    an₁=>2.7, an₂=>1.7, #total factor productivity at max
    γ₁=>0.1, γ₂=>0.1, #relaxation speed towards an₁ and an₂
    δ₁=>0.05, δ₂=>0.05, #capital entropic decay rates  
    s₁=>0.25, s₂=>0.21, #savings rates (of disposable income)
    Cₘ₁=>0.7, Cₘ₂=>0.7, #minimum subsistence consumption
    σ₁=>0.03, σ₂=>0.03, dmₓ=> 100, #climate damages on infrastructure
    re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>20,  #decarb, #carbon intensities
    η=>1, #rate at which decarbonization is initiated
    u=>0.014, #earth system carbon uptake and release
    G₁=>5, G₀=>20, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
    ];
    shock_drift_prob = ibri_shock_drift_model_mod.construct_sde_problem(shock_drift_sde, tspan, shock_drift_u0, shock_drift_p)

    function shock_drift_prob_func(prob, i, repeat)
            #rand_Tg = rand(pa)
            #rand_G_0 = rand(Gms)
            #rand_Tg = pa[i]
            #rand_G_0 = Gms[i]
            pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
            ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, an₁=>2.7, an₂=>1.7,
            δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
            σ₁=>combs1[i][3], σ₂=>combs1[i][3], dmₓ=> 100, re₁=>combs1[i][4], 
            re₂=>combs1[i][4], eb=>0.00004, Tg=>combs1[i][1], 
            η=>1, u=>0.014,  G₁=>5, G₀=>combs1[i][2], Gₘ=>20]
            # push!(rand_values, [rand_Tg, rand_G_0])
            remake(prob, p = pnew) #change Tg and G_0
    end

    # build jump problem
    λ₀=0.05
    k=1
    Gₚ=3
    scaling = 0.5 # average damage = 0.1% of K(100)
    crate1(u, p, t) = ibri_shock_drift_model_mod.expit(u[5], λ₀, k, Gₚ)
    crate2(u, p, t) = ibri_shock_drift_model_mod.expit(u[5], λ₀, k, Gₚ)

    function affect1!(integrator)
        integrator.u[9] -= rand(0:0.01:0.2)*integrator.u[9]
        nothing
    end

    function affect2!(integrator)
        integrator.u[10] -= rand(0:0.01:0.2)*integrator.u[10]
        nothing
    end

    crj1 = ConstantRateJump(crate1, affect1!)
    crj2 = ConstantRateJump(crate2, affect2!)
    shock_drift_jump_prob = JumpProblem(shock_drift_prob, Direct(), crj1, crj2)
    # Unified construct_solve using output_func
    function construct_solve(prob, prob_func, N_samp)
        ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
        return solve(ensemble_prob, trajectories = N_samp, saveat=1, SRIW1())
    end

    # Solve the shock drift jump problem
    return construct_solve(shock_drift_jump_prob, shock_drift_prob_func, N_samp)
end # end of function
end # end of module

# --- Run all models and plot results ---
std_sol = StandardModels.run_std_model()
stochastic_sol = StandardModels.run_stochastic_model()
shock_sol = ShockModels.run_shock_model()
drift_sol = DriftModels.run_std_drift_model()
stochastic_drift_sol = DriftModels.run_stochastic_drift_model()
shock_drift_sol = ShockDriftModels.run_shock_drift_model()

# --- Plot ---
include("../../src/ibri_std_model.jl")
using .ibri_std_model_mod

# define saving locations
save_folder = "./figures"


Tg_v1 = [4.6, 5.4]
G0_v1 = [5.0]
sigma1_v1 = [0.03]
re1_v1 = [0.1]
combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))
N_samp = length(combs1)

# pick 3 categorical colors from batlow
cols = cgrad(:batlow, 4, categorical=true)

labels = [" Favorable", " Drifting", " Critical"]

# wider figure: 1 row × 3 cols
pp = plot(layout=(1,3),
          xlim=(250,750), ylim=(0,170),
          grid=false,
          size=(800,250),
          legend=true,
          title = ["A" "B" "C"], 
          titlelocation = :left)

# --- First panel (left) ---
plot!(pp[1], std_sol[1][G]*100, std_sol[1][Y₁+Y₂],
      color=cols[1], label="")
plot!(pp[1], std_sol[2][G]*100, std_sol[2][Y₁+Y₂],
      color=cols[2], label="")
plot!(pp[1], drift_sol[2][G]*100, drift_sol[2][Y₁+Y₂],
      color=cols[3], label="")
plot!(pp[1], ylabel="World GNI [Trillion \$]")

# --- Second panel (middle) ---
plot!(pp[2], stochastic_sol[1][G]*100, stochastic_sol[1][Y₁+Y₂],
      color=cols[1], label="")
plot!(pp[2], stochastic_sol[2][G]*100, stochastic_sol[2][Y₁+Y₂],
      color=cols[2], label="")
plot!(pp[2], stochastic_drift_sol[2][G]*100, stochastic_drift_sol[2][Y₁+Y₂],
      color=cols[3], label="")
plot!(pp[2], xlabel="Atmospheric Carbon Concentration [ppm]")

# we include the shock model file here to avoid namespace conflicts, note that rerunning can lead to issues --> start new REPL
include("../../src/ibri_shock_model.jl")
using .ibri_shock_model_mod

# --- Third panel (right) ---
plot!(pp[3], shock_sol[1][ibri_shock_model_mod.G]*100,
             shock_sol[1][ibri_shock_model_mod.Y₁+ibri_shock_model_mod.Y₂],
      color=cols[1], label=labels[1])
plot!(pp[3], shock_drift_sol[2][ibri_shock_model_mod.G]*100,
             shock_drift_sol[2][ibri_shock_model_mod.Y₁+ibri_shock_model_mod.Y₂],
      color=cols[3], label=labels[2])
plot!(pp[3], shock_sol[2][ibri_shock_model_mod.G]*100,
             shock_sol[2][ibri_shock_model_mod.Y₁+ibri_shock_model_mod.Y₂],
      color=cols[2], label=labels[3])

# remove redundant ticks
plot!(pp[2], yticks=false)
plot!(pp[3], yticks=false)

# margins
plot!(pp[1], left_margin=5mm, bottom_margin=9mm)

println("Saving...")
savefig(pp, joinpath(save_folder, "figure5.svg"))
savefig(pp, joinpath(save_folder, "figure5.png"))
savefig(pp, joinpath(save_folder, "figure5.pdf"))
println("Done Saving.")