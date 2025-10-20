# © 2025 Max Bechthold, John M. Anderies and the IBRI team

using ModelingToolkit, DifferentialEquations, Plots, Arrow, DataFrames, Statistics, CSV, Plots.PlotMeasures, Colors, ColorSchemes

"""This script plots Figure 6 as found in the publication:
It creates a 2x2 panel plot comparing the resilience index results for the deterministic many worlds with gridded parameters (A),
the stochastic many worlds with shocks (B), the stochastic many worlds with small multiplicative noise (C) and the stochastic many worlds with large multiplicative noise (D).
It also produces a CSV file with resulting overall resilience index values."""

# load models and modules
include("../../src/ibri_std_model.jl")
include("../../src/ibri_analyses.jl")


using .ibri_std_model_mod
using .ibri_analyses_mod

# define data and saving locations
data_folder = "./data"
save_folder = "./figures"
results_folder = "./results"

# load data
println("Loading Data...")
df_det = DataFrame(Arrow.Table(joinpath(data_folder, "deterministic_many_worlds_compare.arrow")))
df_shock = DataFrame(Arrow.Table(joinpath(data_folder, "stochastic_many_worlds_shocks.arrow")))
df_25 = DataFrame(Arrow.Table(joinpath(data_folder, "stochastic_many_worlds_mult_small_sigma.arrow")))
df_75 = DataFrame(Arrow.Table(joinpath(data_folder, "stochastic_many_worlds_mult_large_sigma.arrow")))
println("Done Loading Data.")


# all share the same parameter ranges
Tg_v = collect(4.2:0.025:6.2)
G0_v = collect(4.2:0.025:6.2)
sigma1_v = [0.025, 0.05, 0.075, 0.1]
re1_v = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15]

# define helper function to calculate resilience values (btw, julia docstrings go before the object)
"""
Compute a matrix of mean resilience values for given `Tg` and `G0` parameter ranges.
# Arguments
- `df::DataFrame` : Input DataFrame containing columns `:Tg`, `:G0`, and `:value`.
- `Tg_v::AbstractVector` : Vector of `Tg` values (columns of the output matrix).
- `G0_v::AbstractVector` : Vector of `G0` values (rows of the output matrix).

# Returns
- `res_matrix::Matrix{Float64}` : A matrix of size `(length(G0_v), length(Tg_v))` where each element corresponds
to the mean of `:value` for the `(G0, Tg)` pair. Missing pairs are filled with `NaN`.
"""
function resilience_matrix(df, Tg_v, G0_v)
    df_mean = combine(groupby(df, [:Tg, :G0]), :value => mean => :value_mean)
    res_dict = Dict((row.Tg, row.G0) => row.value_mean for row in eachrow(df_mean))
    res_matrix = [get(res_dict, (g0, tg), NaN) for g0 in G0_v, tg in Tg_v]
    return res_matrix
end

# calculate the matrices
mats = [
    resilience_matrix(df_det,   Tg_v, G0_v),
    resilience_matrix(df_shock, Tg_v, G0_v),
    resilience_matrix(df_25,    Tg_v, G0_v),
    resilience_matrix(df_75,    Tg_v, G0_v)
]

# set titles
panel_titles = ["A", "B", "C", "D"]

# create heatmaps
p1 = heatmap(G0_v*100, Tg_v*100, mats[1], title=panel_titles[1], titlefont=font(18), legend=false, xticks=false, yticks=:auto, tickfont=font(14), c=cgrad(:vik, rev=true), clim=(0,1))
p2 = heatmap(G0_v*100, Tg_v*100, mats[2], title=panel_titles[2], titlefont=font(18), legend=false, xticks=false, yticks=false, tickfont=font(14), c=cgrad(:vik, rev=true), clim=(0,1))
p3 = heatmap(G0_v*100, Tg_v*100, mats[3], title=panel_titles[3], titlefont=font(18), legend=false, xticks=:auto, yticks=:auto, tickfont=font(14), c=cgrad(:vik, rev=true), clim=(0,1))
p4 = heatmap(G0_v*100, Tg_v*100, mats[4], title=panel_titles[4], titlefont=font(18), legend=false, xticks=:auto, yticks=false, tickfont=font(14), c=cgrad(:vik, rev=true), clim=(0,1))

# create colorbar (cb) 
# to bypass the problem that Plots does not support one global colorbar, we create a "colorbar" as a "5th heatmap"
cb_vals = (0:0.001:1) .* ones(1001,1) # need a high resolution here
cb = heatmap(cb_vals, c=cgrad(:vik, rev=true), legend=false, xticks=false, yticks=(1:100:1001, string.(0:0.1:1)), tickfont=font(14), ylabel="Resilience Index IBRI", yguidefontsize=18)

# same problem for shared x and ylabels, so we create these additional spacings
xlabel_spacing = plot(legend=false, framestyle=:none, xaxis=false, yaxis=false, margin=0mm)
annotate!(xlabel_spacing, 0.5, 0.5, text("Climate Threshold G₀ [ppm]", :center, 18)) # add the wanted label text
plot!(xlabel_spacing, bottom_margin=-1mm) # take some unneeded space away

# repeat for y label
ylabel_spacing = plot(legend=false, framestyle=:none, xaxis=false, yaxis=false, margin=0mm)
annotate!(ylabel_spacing, 0.5, 0.5, text("Decarbonizationpoint Tg [ppm]", :center, 18, rotation=90))
plot!(ylabel_spacing, left_margin=-4mm)

# repeat for colorbar
cb_spacing = plot(legend=false, framestyle=:none, xaxis=false, yaxis=false, margin=0mm)

# use layout with the corresponding elements and widths/heigths (in fractions)
l = @layout [
    ylabel_spacing{0.025w}   [grid(2,2); xlabel_spacing{0.035h}]  cb_spacing{0.0001w} cb{0.04w}
]

# finally plot
p = plot(ylabel_spacing, p1, p2, p3, p4, xlabel_spacing, cb_spacing, cb;
         layout = l, size=(1400,950),
         titlelocation=:left)

# save
println("Saving...")
savefig(p, joinpath(save_folder, "figure6.svg"))
savefig(p, joinpath(save_folder, "figure6.png"))
savefig(p, joinpath(save_folder, "figure6.pdf"))
println("Done Saving.")


println("Calculating overall resilience index values...")
# calculate overall resilience index values for both cases and sub cases

function resiliences(df, Tg_v, G0_v)
    df_mean = combine(groupby(df, [:Tg, :G0]), :value => mean => :value_mean)
    res_dict = Dict((row.Tg, row.G0) => row.value_mean for row in eachrow(df_mean))
    total_res = mean(values(res_dict))
    return total_res
end

# List of datasets and labels
datasets = [
    (df_det, "Deterministic"),
    (df_shock, "Stochastic Shocks"),
    (df_25, "Stochastic Small Noise"),
    (df_75, "Stochastic Large Noise")
]

# Compute overall resilience for each dataset
metrics = String[]
resilience_values = Float64[]

for (df, label) in datasets
    push!(metrics, label)
    push!(resilience_values, resiliences(df, Tg_v, G0_v))
end

# Create DataFrame
df_res = DataFrame(
    Metric = metrics,
    Resilience = resilience_values
)

# Save to CSV
CSV.write(joinpath(results_folder, "overall_resilience_values_stochastic.csv"), df_res)
println("Done Calculating.")