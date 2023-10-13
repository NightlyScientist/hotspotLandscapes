using JLD2
using StatsBase
using StaticArrays
import Measurements: ±, value, uncertainty, measurement
using CairoMakie, Colors 
using ColorSchemes
import ColorSchemes: get, hot, amp

include("../../src/tools/parseFiles.jl");
include("../../src/tools/binUtils.jl");
include("../../src/tools/containers.jl")
include("../../src/base/model.jl");
include("../../src/base/dataModels.jl")
include("../../src/base/environment.jl")

using Revise
includet("measures.jl")
includet("figures.jl")

import .GraphModel: gNodes, buildMap, overlap

# check if being ran interactively
scriptMode = (abspath(PROGRAM_FILE) == @__FILE__)

# find data fSets from path
if scriptMode
    addargs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

    global inputs = begin 
        sts = ArgParseSettings()
        addargs!(sts, "path", arg_type=String, required=true)
        addargs!(sts, "smpl", action=:store_true)
        addargs!(sts, "nl", default=5)

        namedtuple(parse_args(sts))
    end
else
    basePath = mkpath(readline("analysis/shortestpaths/output/defaults/basedirectory.txt"))
    sampleSet = (@isdefined pths) ? pths : readlines("analysis/shortestpaths/output/defaults/defaultsample.txt")
    inputs = (nl = 15, smpl = false)
end

# collect simulation ensemble
fSets = append!(map(getCollections, sampleSet)...)

vals = []
for i in eachindex(fSets)
    p = Symbol(first(fSets[i]).parameter)
    push!(vals, getfield.(fSets[i], p))
end

# theme settings
set_theme!(theme_minimal(); figure_padding = 10, fontsize=20)
mytheme = Theme(
    Figure = (backgroundcolor = :transparent),
    Heatmap = (colormap = amp, lowclip = :transparent),
    Scatter = (markersize = 20, colormap = ColorScheme([RGBA(1,1,1, 0.0), RGBA(0.4,0.8,0.0,0.5), RGBA(0,0.5,0,0.5), RGBA(1,1,0,.5)])
    )
);
update_theme!(mytheme)

# containers for processed data
dataModels = Vector{NamedTuple}(undef, length(fSets))

@info "analysis/shortestpaths/output/databases/datamodels_$(pths |> first |> splitdir |> last).jld2"

jldsave("analysis/shortestpaths/output/databases/datamodels_$(pths |> first |> splitdir |> last).jld2"; datamodels = dataModels)

# dataModels = load("analysis/shortestpaths/output/databases/datamodels.jlds2", "datamodels")

# process fSets
let nl = inputs.nl, selector = 15

    for (i, partition) in enumerate(fSets)
        # simulation and prediction input data
        rp = first(partition)
        input = jldopen(rp.path * "/data.jld2", "r") 
        spinput = jldopen(readlink(rp.path * "/spaths.jld2"), "r") 

        if ~isassigned(dataModels, i)
            # construct graph
            xy, cntns = buildMap(rp.dims, gNodes(), Tuple, 
                (n,x,y) -> (sqrt(3)*(x - 0.5*(y%2)), 1 + 1.5*(y-1)))

            # build kdtree from hotspot locations
            kdtree = buildHSKDT(xy, input["htspts"], rp.dims)

            # estimate rand Selection choices
            bckgrndEstmt = backgroundEstimate(spinput, input, xy, rp.dims, rp.density, rp.radius; nl = nl)

            # calculate shortest path tracks and visited locations
            trails, roots, sphs = spTracks(spinput, xy, rp.dims, kdtree, rp.radius; nl = nl)

            # predict hs visited by SPs
            acc_hs, freq_hs = subsetHotspots(kdtree, xy, rp.radius, sphs, input)

            # predict ancestors
            acc_anc, freq_anc, acc_anc_auc = subsetAncestors(input, rp.dims, roots)

            dataModels[i] = (
                config = rp,
                parameter = measurement(mean_and_std(vals[i])...),
                acc_ancestors = acc_anc,
                acc_hotspots = acc_hs,
                freq_hotspots = freq_hs,
                freq_ancestors = freq_anc,
                acc_hotspots_bg = bckgrndEstmt,
                acc_ancestors_auc = acc_anc_auc,
                # freq_survival_anc = extra,
                trails = trails,
                roots = roots
            )
        else
            println("already calculated, just retrieving.")
        end

        # genetic lineage heat map
        lhmap, _ = lineageTracks(rp.dims, input)

        # overlay snapshots 
        with_theme(mytheme) do 

            # fig = trailsImg(input["env"], lhmap, trials, dataModels[i])
            fig = trailsImg(input["env"], lhmap, dataModels[i].trails, dataModels[i])

            survivalFreqImg!(fig, dataModels[i])

            fname = splitdir(rp.path) |> last

            save("analysis/shortestpaths/output/figures/optpaths/left_vf_w_sp_$(fname).png", fig)
            if selector == i || selector == 0
                scriptMode || display(fig)
            end
        end

        close(spinput)
        close(input)
    end
end


# agreement accuracy figure
begin 
    fig = accImg(dataModels)
    # save("analysis/shortestpaths/output/predicitonAccuracy.png", fig)
    scriptMode || display(fig)
end

dm_1 = load("/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/analysis/shortestpaths/output/databases/datamodels_dims:1000.0,1000.0_hsRegion:None_I:xxx_D:0.1_R:10_trials:1200_nCollections:1_intervals:0.0,0.05,1.0_landscape:varied.jld2", "datamodels")

dm_2 = load("/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/analysis/shortestpaths/output/databases/datamodels_dims:1000.0,1000.0_hsRegion:None_I:xxx_D:0.1_R:10_trials:1200_nCollections:1_intervals:0.0,1.0,8.0_landscape:varied.jld2", "datamodels")

let 
    fig = Figure(resolution=(800,750))
    ax = Axis(fig[1, 1], 
        ylabel = "Overlap Fraction", xlabel = L"\phi",
        xlabelsize = 40, ylabelsize = 28,
        yticklabelsize = 25, xticklabelsize = 25
    )

    colors = colorscheme_alpha(roma, 0.5)
    for (i, dms) in enumerate([dm_1, dm_2])
        # hotspots -> agreement with geoemtric paths
        x, xerr = getfield.(dms, :parameter) |> unzip
        y, yerr = getfield.(dms, :acc_hotspots) |> unzip

        # draw accuracy graph for hotspots
        scatter!(ax, x, y, markersize = 35, color=get(colors, 1), 
            label = "Hotspots (M)"
        )
        errorbars!(ax, x, y, xerr, color=get(colors, 1), direction=:x)
        errorbars!(ax, x, y, yerr, color=get(colors, 1))

        # estimates of background for hotspot measure
        bg, bgerr = getfield.(dms, :acc_hotspots_bg) |> unzip

        if i == 2
            lines!(ax, x, fill(value(first(bg)), length(x)),
                color=:red, 
            )
            text!(ax, "Random Hotspot Selection", position = (2.5, 0.34), align = (:left, :center),rotation = 0.0, fontsize = 23)
        end

        # ancestors -> agreement with geoemtric paths
        y, yerr = getfield.(dms, :acc_ancestors) |> unzip

        # draw accuracy graph for ancestors
        scatter!(ax, x[2:end], y[2:end], color=get(colors, 0),
            label = "Ancestors (K)", markersize = 35
        )
        errorbars!(ax, x[2:end], y[2:end], xerr[2:end], color=get(colors, 0), direction=:x)
        errorbars!(ax, x[2:end], y[2:end], yerr[2:end], color=get(colors, 0))
    end
    
    ylims!(ax, (0.2, 1.05))
    axislegend(ax, position=:rb, labelsize=25, merge = true)
    scriptMode || display(fig)
    save("analysis/shortestpaths/output/figures/hotspot_ancestor_SP_overlap.png", fig)
end

let
    fig = Figure()
    ax = Axis(fig[1, 1], backgroundcolor = :transparent)
    rowsize!(fig.layout, 1, Relative(9 / 10))
    # heatmap!(ax, 
    #     reshape(repeat(dataModels[27].freq_ancestors,2), (length(dataModels[27].freq_ancestors), 2)), lowclip=:white
    # )

    @info dataModels[27].config.intensity
    lines!(dataModels[27].freq_ancestors)
    lines!(dataModels[27].roots * maximum(dataModels[27].freq_ancestors))
    
    @info sum(dataModels[27].roots .* dataModels[27].freq_ancestors) / sum(dataModels[27].freq_ancestors)
    # Colorbar(fig[3, 1]; label = "Survival Probability", 
    #     ticksize=20, tickalign=1, ticklabelsize = 22,
    #     colormap = amp,
    #     vertical = false, flipaxis = false, labelsize = 35
    # )

    # hidedecorations!(ax; grid = true)
    xlims!(ax, (1,length(dataModels[27].freq_ancestors)))
    # rowgap!(fig.layout, 3)
    scriptMode || display(fig)
end