# © 2025 Max Bechthold, John M. Anderies and the IBRI team

module ibri_analyses_mod # creates the module

using ModelingToolkit, DifferentialEquations # relevant imports

# includes the file from the model ibri for the namespace
include("ibri_std_model.jl")


"""Calculates the resilience index for a solution sol and a given resilience condition.

sol: the solution of the ODE or SDE problem, where the output function returns G and Y/P for both regions at the last time step.
res_condition: array, first entry gives ppm threshold, second one gives gdp per capita

Btw. Julia Docstrings go before the function definition, not inside like in Python."""
function safe_and_just_resilience_index(sol, res_condition)
    R_i=[]
    for i in eachindex(sol)
        if sol[i][1] < res_condition[1] && sol[i][2] > res_condition[2] && sol[i][3] > res_condition[2]
            push!(R_i, 1)
        else
            push!(R_i, 0)
        end
    end
    return R_i
end


"""Solves an ODE set for 6 parameters and their combinations. Post-processes.
Here, sigma1 = sigma2 and re1 = re2.

prob: the ODE prob created with ModelingToolkit
combs: the combinations of relevante sweeping parameters (σ, re, Tg and G0)
res_condition: array, first entry gives ppm threshold, second one gives gdp per capita"""
function many_worlds_analysis_same_regions_safe_and_just(prob, combs, res_condition=[3.5, 4])

    N_samp = length(combs)
    println("$N_samp runs.")
    function many_worlds_prob_func(prob, i, repeat)
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs[i][3], σ₂=>combs[i][3], dmₓ=> 100, re₁=>combs[i][4], 
        re₂=>combs[i][4], eb=>0.00004, Tg=>combs[i][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs[i][2], Gₘ=>20]
        remake(prob, p = pnew) 
    end # end of function
    output_func(sol, i) = ([sol[end][5], sol[Y₁/P₁][end], sol[Y₂/P₂][end]], false) # we do this to only save the last value of the ppm and reduce solution size to prevent the OOM error
    ensemble_prob = EnsembleProblem(prob, prob_func = many_worlds_prob_func, output_func=output_func)
    sol = solve(ensemble_prob, EnsembleThreads(), trajectories = N_samp)    
    dict = Dict(combs .=> safe_and_just_resilience_index(sol, res_condition))
    return dict
end # end of function


"""Analysis where G0 is somehow drawn from a distribution.
Must be passed as p1, i.e. an Array. The others are sweeped.

prob: the ODE prob
p1: array of G0 values drawn from a distribution
p2, p3, p4: arrays of the other sweeping parameters (Tg, sigma1, re1)
res_condition: array, first entry gives ppm threshold, second one gives gdp per capita"""
function threshold_stochastic_parameters_analysis_safe_and_just(prob, p1, p2, p3, p4, res_condition)
    N_ens = length(p1)
    combs = collect(Iterators.product(p2, p3, p4))
    N_traj = length(combs) * N_ens
    N_combs = length(combs)

    ensemble_indices = 1:N_combs

    ensemble_limits = [0]
    for i in 1:N_traj
        if i % N_ens == 0
            push!(ensemble_limits, i)
        end
    end

    function ensemble_prob_func(prob, i, repeat)
        """Had to build this ugly for/if control sequence since prob_func essentially
        works like a for loop but is a function of which i cannot alter the args,
        such that it was impossible to build a normal count to know in which batch one
        is. The way it works now: i gets counted for every run of the solver so i in
        [1, N_traj] and i compare it the to the bounds of the ensemble_limits."""
        count = 1
        for j in 1:(length(ensemble_limits)-1)
          if i > ensemble_limits[j] && i <= ensemble_limits[j+1] 
            count = ensemble_indices[j]
          elseif i == last(ensemble_limits)
            count = last(ensemble_indices)
          end
        end
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs[count][2], σ₂=>combs[count][2], dmₓ=> 100, re₁=>combs[count][3], 
        re₂=>combs[count][3], eb=>0.00004, Tg=>combs[count][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>p1[i-((count-1)*N_ens)], Gₘ=>20]
        # println(combs[count][1])
        remake(prob, p = pnew)
    end
    output_func(sol, i) = ([sol[end][5], sol[Y₁/P₁][end], sol[Y₂/P₂][end]], false)
    ensemble_prob = EnsembleProblem(prob, prob_func = ensemble_prob_func, output_func=output_func)
    sol = solve(ensemble_prob, EnsembleDistributed(), trajectories = N_traj)    
    sol2 = zeros(N_combs)
    for i in 1:N_traj
        for j in 1:(length(ensemble_limits)-1)
            if i > ensemble_limits[j] && i <= ensemble_limits[j+1] 
                if sol[i][1] < res_condition[1] && sol[i][2] > res_condition[2] && sol[i][3] > res_condition[2]
                    sol2[j] =  sol2[j] + 1
                else
                    sol2[j] = sol2[j] + 0 
                end
            end
        end
    end
    sol3 = sol2./N_ens
    new_combs=[]
    for i in eachindex(combs)
        push!(new_combs, combs[i])
    end
    dict = Dict(new_combs .=> sol3)
    return dict
end 


"""Do the resilience index analysis for a SDE problem, but here, sigma1 = sigma2 and re1 = re2.

prob: the SDE prob
N_ens: number of ensemble members per combination of parameters
p1, p2, p3, p4: arrays of the sweeping parameters (G0, Tg, sigma1, re1)
res_condition: array, first entry gives ppm threshold, second one gives gdp per capita"""
function ensemble_sde_analysis_same_regions_safe_and_just(prob, N_ens, p1, p2, p3, p4, res_condition)
    combs = collect(Iterators.product(p1, p2, p3, p4))
    N_traj = length(combs) * N_ens
    N_combs = length(combs)

    ensemble_indices = 1:N_combs

    ensemble_limits = [0]
    for i in 1:N_traj
        if i % N_ens == 0
            push!(ensemble_limits, i)
        end
    end

    println("$N_traj runs.")

    function ensemble_prob_func(prob, i, repeat)
        count = 1
        for j in 1:(length(ensemble_limits)-1)
            if i > ensemble_limits[j] && i <= ensemble_limits[j+1] 
              count = ensemble_indices[j]
            elseif i == last(ensemble_limits)
              count = last(ensemble_indices)
            end
          end
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs[count][3], σ₂=>combs[count][3], dmₓ=> 100, re₁=>combs[count][4], 
        re₂=>combs[count][4], eb=>0.00004, Tg=>combs[count][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs[count][2], Gₘ=>20]
        remake(prob, p = pnew)
    end
    output_func(sol, i) = ([sol[end][5], sol[Y₁/P₁][end], sol[Y₂/P₂][end]], false) 
    ensemble_prob = EnsembleProblem(prob, prob_func = ensemble_prob_func, output_func=output_func)
    println("Start solving...")
    sol = solve(ensemble_prob, EnsembleThreads(), trajectories = N_traj)
    println("Solving done!")
    sol2 = zeros(N_combs)
    for i in 1:N_traj
        for j in 1:(length(ensemble_limits)-1)
            if i > ensemble_limits[j] && i <= ensemble_limits[j+1]
                if sol[i][1] < res_condition[1] && sol[i][2] > res_condition[2] && sol[i][3] > res_condition[2]
                    sol2[j] =  sol2[j] + 1
                else
                    sol2[j] = sol2[j] + 0 
                end
            end
        end
    end
    sol3 = sol2./N_ens
    new_combs=[]
    for i in eachindex(combs)
        push!(new_combs, combs[i])
    end
    dict = Dict(new_combs .=> sol3)
    return dict
end 


"""Do the resilience index analysis for a jump SDE problem.

prob: the jump prob
N_ens: number of ensemble members per combination of parameters
p1, p2, p3, p4: arrays of the sweeping parameters (G0, Tg, sigma1, re1)
res_condition: array, first entry gives ppm threshold, second one gives gdp per capita"""
function ensemble_jump_analysis_same_regions_safe_and_just(prob, N_ens, p1, p2, p3, p4, res_condition)
    combs = collect(Iterators.product(p1, p2, p3, p4))
    N_traj = length(combs) * N_ens
    N_combs = length(combs)

    ensemble_indices = 1:N_combs

    ensemble_limits = [0]
    for i in 1:N_traj
        if i % N_ens == 0
            push!(ensemble_limits, i)
        end
    end

    println("$N_traj runs.")

    function ensemble_prob_func(prob, i, repeat)
        count = 1
        for j in 1:(length(ensemble_limits)-1)
            if i > ensemble_limits[j] && i <= ensemble_limits[j+1] 
              count = ensemble_indices[j]
            elseif i == last(ensemble_limits)
              count = last(ensemble_indices)
            end
          end
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, an₁=>2.7, an₂=>1.7,
        δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs1[i][3], σ₂=>combs1[i][3], dmₓ=> 100, re₁=>combs1[i][4], 
        re₂=>combs1[i][4], eb=>0.00004, Tg=>combs1[i][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2], Gₘ=>20]
        remake(prob, p = pnew)
    end
    output_func(sol, i) = ([sol[end][5], sol[Y₁/P₁][end], sol[Y₂/P₂][end]], false) # we do this to only save the last value of the ppm and reduce solution size to prevent the OOM error
    ensemble_prob = EnsembleProblem(prob, prob_func = ensemble_prob_func, output_func=output_func)
    println("Start solving...")
    sol = solve(ensemble_prob,  SRIW1(), EnsembleThreads(), trajectories = N_traj)
    println("Solving done!")
    sol2 = zeros(N_combs)
    for i in 1:N_traj
        for j in 1:(length(ensemble_limits)-1)
            if i > ensemble_limits[j] && i <= ensemble_limits[j+1]
                if sol[i][1] < res_condition[1] && sol[i][2] > res_condition[2] && sol[i][3] > res_condition[2]
                    sol2[j] =  sol2[j] + 1
                else
                    sol2[j] = sol2[j] + 0 
                end
            end
        end
    end
    sol3 = sol2./N_ens
    new_combs=[]
    for i in eachindex(combs)
        push!(new_combs, combs[i])
    end
    dict = Dict(new_combs .=> sol3)
    return dict
end 

"""Util Function to calculate excess CO2 above pre-industrial.

concentration: CO2 concentration in 100 ppm"""
function excess_CO2(concentration)

    base_concentration = 2.8
    return concentration - base_concentration
end


"""Saves a config.

PATH = path where to save
prob = ODE or SDEProblem (this contains tspan etc. as well, but in a strange format)
tspan = timeframe
u0 = initial conditions
p = parameters
ranges = parameter ranges that one looped
further_params = anything else you want to save (neednt be, i.e. nothing)"""
function save_config(PATH, prob, tspan, u0, p, ranges, further_params)

    jldsave("$PATH/config.jld2"; prob, tspan, u0, p, ranges, further_params)
end

end # end of module