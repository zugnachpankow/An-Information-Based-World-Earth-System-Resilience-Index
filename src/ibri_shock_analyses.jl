# © 2025 Max Bechthold, John M. Anderies and the IBRI team

module ibri_shock_analyses_mod
# creates the moduls

# relevant imports
using ModelingToolkit, DifferentialEquations

# includes the model from the model ibri
include("ibri_shock_model.jl")


#-----definition of relevant functions-----


function ensemble_jump_analysis_same_regions_safe_and_just(prob, N_ens, p1, p2, p3, p4, res_condition)
    """
    Do the resilience index analysis for a jump SDE problem, but
    here, sigma1 = sigma2 and re1 = re2.
    lambda, k, Gp are the parameters for the expit
    """
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
        """Had to build this ugly for/if control sequence since prob_func essentially
        works like a for loop but is a function of which i cannot alter the args,
        such that it was impossible to build a normal count to know in which batch one
        is. The way it works now: i gets counted for every run of the solver so i in
        [1, N_traj] and i compare it the to the bounds of the ensemble_limits."""
        # println("$i of $N_traj...")
        
        count = 1
        # print(count)
        for j in 1:(length(ensemble_limits)-1)
            # println(k)
            if i > ensemble_limits[j] && i <= ensemble_limits[j+1] 
              # println(i)
              #println(ensemble_limits[k], ensemble_limits[k+1], ensemble_indices[k])
              count = ensemble_indices[j]
            elseif i == last(ensemble_limits)
              count = last(ensemble_indices)
            end
            # println(i)
          end
        # println(count)
        # println(ensemble_limits[repeat])
        # if i > ensemble_limits[repeat]
        #   repeat += 1
        # end
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, an₁=>2.7, an₂=>1.7,
        δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs[count][3], σ₂=>combs[count][3], dmₓ=> 100, re₁=>combs[count][4], 
        re₂=>combs[count][4], eb=>0.00004, Tg=>combs[count][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs[count][2], Gₘ=>20]
        # println(combs[count][1])
        remake(prob, p = pnew)
    end
    output_func(sol, i) = ([sol[end][5], sol[Y₁/P₁][end], sol[Y₂/P₂][end]], false) # we do this to only save the last value of the ppm and reduce solution size to prevent the OOM error
    ensemble_prob = EnsembleProblem(prob, prob_func = ensemble_prob_func, output_func=output_func)
    println("Start solving...")
    sol = solve(ensemble_prob,  Tsit5(), EnsembleThreads(), trajectories = N_traj)
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

end # end of module