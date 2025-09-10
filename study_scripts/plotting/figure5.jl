using ModelingToolkit, DifferentialEquations, Plots, JLD2, DataFrames, Arrow, Random, Colors, Plots.PlotMeasures

include("../../src/ibri_std_model.jl")
using .ibri_std_model_mod

# define data and saving locations
data_folder = "./data"
save_folder = "./figures"

@load "/home/maxbecht/WES_Resilience/code/paper_studies/Results3/std_sol.jld2" std_sol
@load "/home/maxbecht/WES_Resilience/code/paper_studies/Results3/drift_sol.jld2" drift_sol
@load "/home/maxbecht/WES_Resilience/code/paper_studies/Results3/stochastic_sol.jld2" stochastic_sol
@load "/home/maxbecht/WES_Resilience/code/paper_studies/Results3/stochastic_drift_sol.jld2" stochastic_drift_sol
@load "/home/maxbecht/WES_Resilience/code/paper_studies/Results3/shock_sol.jld2" shock_sol
@load "/home/maxbecht/WES_Resilience/code/paper_studies/Results3/shock_drift_sol.jld2" shock_drift_sol

# --- plot ---

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

# we include the shock model file here to avoid namespace conflicts
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

savefig(pp, "/mnt/c/Users/maxbecht/Desktop/WES_paper/figures_paper/stochastic_models.png")
savefig(pp, "/mnt/c/Users/maxbecht/Desktop/WES_paper/figures_paper/stochastic_models.svg")
savefig(pp, "/mnt/c/Users/maxbecht/Desktop/WES_paper/figures_paper/stochastic_models.pdf")