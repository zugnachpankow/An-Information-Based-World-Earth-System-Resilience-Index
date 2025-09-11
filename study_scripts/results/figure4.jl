# © 2025 Max Bechthold, John M. Anderies and the IBRI team

using ModelingToolkit, DifferentialEquations, Arrow, DataFrames, CSV, Statistics, Plots.PlotMeasures, Colors, ColorSchemes

"""This script plots Figure 4 as found in the publication:
It creates a plot showing the IBRI resilience index as a function of decarbonization point Tg and climate threshold G0 (A)
and as a function of decarbonization point Tg and decarbonization rate re1 (B) with a reailistically distributed threshold.
It also produces a CSV file with resulting overall resilience index values."""

include("../../src/ibri_std_model.jl")
include("../../src/ibri_analyses.jl")

using .ibri_std_model_mod
using .ibri_analyses_mod

# define data and saving locations
data_folder = "./data"
save_folder = "./figures"
results_folder = "./results"

println("Loading Data...")
df_g = DataFrame(Arrow.Table(joinpath(data_folder, "deterministic_many_worlds_gridded.arrow")))
df_s = DataFrame(Arrow.Table(joinpath(data_folder, "deterministic_many_worlds_threshold_distribution.arrow")))
println("Done Loading Data.")

# full info set, is same for both
Tg_v = collect(4.2:0.01:6.2)
G0_v = collect(4.2:0.01:6.2)
sigma1_v = collect(0.01:0.01:0.1)
re1_v = collect(0.01:0.01:0.15)
re1_s = collect(0.0:0.005:0.1) # this one differs for the threshold distribution case

# compute mean values for each parameter combination
# and then reshape into matrix for plotting
df_mean_g = combine(groupby(df_g, [:Tg, :G0]), :value => mean => :value_mean)
res_dict_g = Dict((row.Tg, row.G0) => row.value_mean for row in eachrow(df_mean_g))
res_matrix_g = [get(res_dict_g, (g0, tg), NaN) for g0 in G0_v, tg in Tg_v]

df_mean_s = combine(groupby(df_s, [:Tg, :re1]), :value => mean => :value_mean)
res_dict_s = Dict((row.Tg, row.re1) => row.value_mean for row in eachrow(df_mean_s))
res_matrix_s = [get(res_dict_s, (tg, re), NaN) for tg in Tg_v, re in re1_s]

println("Plotting...")

# set titles
panel_titles = ["A", "B"]

# heatmap, moving from fractionals to percentages
p1 = heatmap(G0_v*100, Tg_v*100, res_matrix_g.* 100,  # note the transpose!
    c = cgrad(:vik, rev=true),
    legend = false,
    title = panel_titles[1],
    titlefont=font(18),
    xticks=:auto, 
    yticks=:auto, 
    tickfont=font(14), 
    clim = (0,100),
    xlabel = "Climate Threshold G₀ [ppm]",
    xguidefontsize=18,
    ylabel = "Decarbonizationpoint Tg [ppm]",
    yguidefontsize=18,
    xlims = (minimum(G0_v)*100, maximum(G0_v)*100),
    ylims = (minimum(Tg_v)*100, maximum(Tg_v)*100)
)


plot!([450, 550, 550, 450, 450],
      [450, 450, 550, 550, 450],
      color=:white, lw=2, linestyle=:solid, label="R1")

plot!([440, 520, 520, 440, 440],
      [470, 470, 520, 520, 470],
      color=:white, lw=2, linestyle=:dash, label="R2")

# i need to bypass the standard legend command as it will trigger a colorbar
legend_x = minimum(G0_v)*100 + 20  # 20 units inset from the left
legend_y = maximum(Tg_v)*100 - 10  # 10 units inset from the top

annotate!(legend_x + 25 + 4, legend_y - 4, text("Solid: Resilience Set R1", 12, :white))
annotate!(legend_x + 29.5 + 4, legend_y - 16, text("Dashed: Resilience Set R2", 12, :white))

p2 = heatmap(re1_s*100, Tg_v*100, res_matrix_s .* 100,
    c = cgrad(:vik, rev=true),
    legend = false,
    title = panel_titles[2],
    titlefont=font(18),
    xticks=:auto, 
    yticks=false, 
    tickfont=font(14),
    clim = (0,100),
    xlabel = "Decarbonization Efforts per Region [%/yr]",
    xguidefontsize=18,
    ylabel = "",
    xlims = (0.0, 10),
    ylims = (minimum(Tg_v)*100, maximum(Tg_v)*100)
)

# create colorbar (cb) 
# to bypass the problem that Plots does not support one global colorbar, we create a "colorbar" as a "5th heatmap"
cb_vals = (0:0.1:100) .* ones(1001,1) # need a high resolution here
cb = heatmap(cb_vals, c=cgrad(:vik, rev=true), legend=false, xticks=false, yticks=(1:100:1001, string.(0:10:100)), tickfont=font(14), ylabel="Resilience Index IBRI [%]", yguidefontsize=18)

# use layout with the corresponding elements and widths/heigths (in fractions)
l = @layout [
    grid(1, 2) cb{0.04w}
]

# finally plot
p = plot(p1, p2, cb;
         layout = l, size=(1400,600),
         titlelocation=:left, left_margin=7mm, bottom_margin=9mm)

println("Done Plotting...")


println("Saving...")
savefig(p, joinpath(save_folder, "figure4.svg"))
savefig(p, joinpath(save_folder, "figure4.png"))
savefig(p, joinpath(save_folder, "figure4.pdf"))
println("Done Saving.")

println("Calculating overall resilience index values...")
# calculate overall resilience index values for both cases and sub cases

# gridded
total_res_g = mean(values(res_dict_g))

# set R1
Tg_r1 = collect(4.5:0.01:5.5)
G0_r1 = collect(4.5:0.01:5.5)
res_r1 = mean([get(res_dict_g, (tg, g0), NaN) for g0 in G0_r1, tg in Tg_r1])

# set R2
Tg_r2 = collect(4.7:0.01:5.2)
G0_r2 = collect(4.4:0.01:5.2)
res_r2 = mean([get(res_dict_g, (tg, g0), NaN) for g0 in G0_r2, tg in Tg_r2])

# threshold distribution
total_res_s = mean(values(res_dict_s))

# 5% decarbonization rate
Tg_r3 = collect(4.2:0.01:6.2)
re_r3 = [0.05]
res_r3 = mean([get(res_dict_s, (tg, re), NaN) for re in re_r3, tg in Tg_r3])

# save to CSV
df_res = DataFrame(
    Metric = ["Total Gridded", "R1", "R2", "Total Threshold Distribution", "5% Decarb Rate"],
    Resilience = [total_res_g, res_r1, res_r2, total_res_s, res_r3]
)

CSV.write(joinpath(results_folder, "overall_resilience_values_deterministic.csv"), df_res)
println("Done Calculating.")