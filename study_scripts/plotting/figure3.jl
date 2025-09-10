using ModelingToolkit, DifferentialEquations, Plots, DataFrames, Arrow, Random, Colors

"""This script plots Figure 3 as found in the publication:
It creates a plot showing the Safe Operating Space (SOS) and Safe and Just Operating Space (SAJOS)"""

include("../../src/ibri_std_model.jl")
include("../../src/ibri_analyses.jl")


using .ibri_std_model_mod
using .ibri_analyses_mod

# define data and saving locations
data_folder = "./data"
save_folder = "./figures"

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
re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.2,  #decarb, #carbon intensities
η=>1, #rate at which decarbonization is initiated
u=>0.014, α=>0.1, #earth system carbon uptake and release
G₁=>5, G₀=>4.6, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
];

# use helper function to construct the ODE problem
prob = ibri_std_model_mod.construct_ode_problem(ode, tspan, u0, p)

# construct arrays for changing params
Tg_v1 = [3.8, 4.2, 4.6, 5.0, 5.4, 5.8, 6.2, 6.6, 7.0] # decarbonization points
G0_v1 = [10] # threshold is so high it is effectively not there
sigma1_v1 = [0.025] # low climate damages
re1_v1 = [0.1] # rather high decarbonization rate
combs1 = unique(collect(Iterators.product(Tg_v1, G0_v1, sigma1_v1, re1_v1)))

# construct different set with higher climate damages
Tg_v2 = [3.8, 4.2, 4.6, 5.0, 5.4, 5.8, 6.2, 6.6, 7.0]
G0_v2 = [10] 
sigma1_v2 = [0.075] #set differs in here
re1_v2 = [0.1]
combs2 = unique(collect(Iterators.product(Tg_v2, G0_v2, sigma1_v2, re1_v2)))

# now we repeat for a low (i.e. "existing") climate threshold
Tg_v3 = [3.8, 4.2, 4.6, 5.0, 5.4, 5.8, 6.2, 6.6, 7.0]
G0_v3 = [5.0] # 
sigma1_v3 = [0.025]
re1_v3 = [0.1]
combs3 = unique(collect(Iterators.product(Tg_v3, G0_v3, sigma1_v3, re1_v3)))

Tg_v4 = [3.8, 4.2, 4.6, 5.0, 5.4, 5.8, 6.2, 6.6, 7.0]
G0_v4 = [5.0]
sigma1_v4 = [0.075]
re1_v4 = [0.1]
combs4 = unique(collect(Iterators.product(Tg_v4, G0_v4, sigma1_v4, re1_v4)))
println("Model setup done!")


println("Start runs...")
# we build multiple prob_funcs to avoid namespace issues
function prob_func1(prob, i, repeat)
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs1[i][3], σ₂=>combs1[i][3], dmₓ=> 100, re₁=>combs1[i][4], 
        re₂=>combs1[i][4], eb=>0.00004, Tg=>combs1[i][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs1[i][2], Gₘ=>20]
        remake(prob, p = pnew) #change Tg and G_0
end

function prob_func2(prob, i, repeat)
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs2[i][3], σ₂=>combs2[i][3], dmₓ=> 100, re₁=>combs2[i][4], 
        re₂=>combs2[i][4], eb=>0.00004, Tg=>combs2[i][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs2[i][2], Gₘ=>20]
        remake(prob, p = pnew) #change Tg and G_0
end

function prob_func3(prob, i, repeat)
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs3[i][3], σ₂=>combs3[i][3], dmₓ=> 100, re₁=>combs3[i][4], 
        re₂=>combs3[i][4], eb=>0.00004, Tg=>combs3[i][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs3[i][2], Gₘ=>20]
        remake(prob, p = pnew) #change Tg and G_0
end

function prob_func4(prob, i, repeat)
        pnew = [r₁=> 0.038, Kp₁=>1.5, r₂=>0.042, Kp₂=>9.7, ml₁₂=>10, 
        ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, r₁₂=>0, r₂₁=>0, α₁=>0.5, α₂=>0.5, a₁=>2.7, 
        a₂=>1.7, δ₁=>0.05, δ₂=>0.05, s₁=>0.25, s₂=>0.21, Cₘ₁=>0.7, Cₘ₂=>0.7, 
        σ₁=>combs4[i][3], σ₂=>combs4[i][3], dmₓ=> 100, re₁=>combs4[i][4], 
        re₂=>combs4[i][4], eb=>0.00004, Tg=>combs4[i][1], 
        η=>1, u=>0.014, α=>0.1, G₁=>5, G₀=>combs4[i][2], Gₘ=>20]
        remake(prob, p = pnew) #change Tg and G_0
end

N_samp1 = length(combs1)
N_samp2 = length(combs2)
N_samp3 = length(combs3)
N_samp4 = length(combs4)
    
# put the solver into a function, since according to documentation thats faster
function construct_solve(prob, prob_func, N_samp)
        ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
        return solve(ensemble_prob, trajectories = N_samp, saveat=1)
end

# no we solve the 4 different sets
sol1 = construct_solve(prob, prob_func1, N_samp1)
sol2 = construct_solve(prob, prob_func2, N_samp2)
sol3 = construct_solve(prob, prob_func3, N_samp3)
sol4 = construct_solve(prob, prob_func4, N_samp4)
println("Runs done!")


# plot SOS, SAJOS and certainly unsafe

# SOS
x_rect1 = [200, 350, 350, 200]
y_rect1 = [0, 0, 170, 170] 

# SAJOS
x_rect2 = [200, 350, 350, 200]
y_rect2 = [100, 100, 170, 170] 

cmap1 = cgrad(:Dark2_4, N_samp1*4, categorical=true, rev=false)

pp = plot(grid=false,
          xlabel="Atmospheric Carbon Concentration [ppm]",
          ylabel="World Gross National Income [Trillion \$]",
          xlim=(200, 800), ylim=(0, 170),
          legend=false)

# SOS
plot!(pp, x_rect1, y_rect1,
      seriestype=:shape,
      c=:lightblue,
      alpha=0.2,
      label="")

# SAJOS
plot!(pp, x_rect2, y_rect2,
      seriestype=:shape,
      c=:lightblue,
      alpha=0.3,
      label="")

# USO
# number of gradient steps
n_steps = 100
x_min, x_max = 500, 800
y_min, y_max = 0, 170

# we try do indicate a region that gets more and more dangerous with increasing G
for i in 1:n_steps
    # interpolation factor from 0 to 1
    t = (i - 1) / (n_steps - 1)

    # color gradient from light red to dark red (or alpha 0 → 0.4)
    color = RGBA(1.0, 0.0, 0.0, 0.0 + t * 0.4)  # keep red, vary alpha

    # strip coordinates
    x1 = x_min + (i - 1) * (x_max - x_min) / n_steps
    x2 = x_min + i * (x_max - x_min) / n_steps
    xs = [x1, x2, x2, x1]
    ys = [y_min, y_min, y_max, y_max]

    plot!(pp, xs, ys,
          seriestype=:shape,
          c=color,
          linealpha=0,
          label="")
end

# add vertical dashed line
plot!(pp, [350, 350], [0, 170],
      linestyle=:dash,
      color=:darkblue,
      linewidth=2,
      label="")

# add horizontal dashed line
plot!(pp, [200, 350], [100, 100],
      linestyle=:dash,
      color=:darkblue,
      linewidth=2,
      label="")

# plot sol1
for i in 1:N_samp1
    plot!(pp, sol1[i][G]*100, sol1[i][Y₁] .+ sol1[i][Y₂],
          alpha=0.9, c=cmap1[i])
end

# plot sol2
for i in 1:N_samp2
    plot!(pp, sol2[i][G]*100, sol2[i][Y₁] .+ sol2[i][Y₂],
          alpha=0.9, c=cmap1[i+N_samp1])
end

# plot sol3
for i in 1:N_samp3
    plot!(pp, sol3[i][G]*100, sol3[i][Y₁] .+ sol3[i][Y₂],
          alpha=0.9, c=cmap1[i+2*N_samp1])
end

# plot sol4
for i in 1:N_samp4
    plot!(pp, sol4[i][G]*100, sol4[i][Y₁] .+ sol4[i][Y₂],
          alpha=0.9, c=cmap1[i+3*N_samp1])
end

# Add text annotation in the middle of the rectangle
annotate!(pp, 275, 135, text("Safe\nAnd\nJust\nOperating\nSpace", 12, :darkblue, :center))
annotate!(pp, 275, 50, text("Safe\nOperating\nSpace", 12, :darkblue, :center))
annotate!(pp, 650, 85, text("Certainly\nUnsafe\nOperating\nSpace", 12, :darkred, :center))

display(pp)
println("Done Plotting...")

println("Saving...")
savefig(pp, joinpath(save_folder, "figure3.svg"))
savefig(pp, joinpath(save_folder, "figure3.png"))
savefig(pp, joinpath(save_folder, "figure3.pdf"))
println("Done Saving.")
