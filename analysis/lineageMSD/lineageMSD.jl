import Base.Threads: @spawn, fetch
import JLD2: jldopen, jldsave
import FileIO: load
import StatsBase: mean, kurtosis, std, median, mode
import Colors: RGBA
import ColorSchemes: get, amp, roma, hot, ColorScheme
import Measurements: value, uncertainty, ±, measurement
import EasyFit.LsqFit: curve_fit
using CairoMakie
using DataStructures
using LaTeXStrings

include("../../src/tools/parseFiles.jl");
include("../../src/tools/binUtils.jl");
include("../../src/base/model.jl");

import .GraphModel: gNodes, buildMap

function colorscheme_alpha(
    cscheme::ColorScheme,
    alpha::T = 0.5;
    ncolors = 12,
) where {T<:Real}
    return ColorScheme([
        RGBA(get(cscheme, k), alpha) for k in range(0, 1, length = ncolors)
    ])
end

# using P.B.C get the arc distance between two points
function periodicBoundaries(halfwidth::Float64, width::Float64)
    f = let h=halfwidth, w=width
        x -> (x > h) ? w - x : x
    end
    return f
end;

function passFltr(y , x = collect(0:length(y)-1); uplim = 1000, offset = 1)
    ps = (x .!== NaN) .& (y .!== NaN) .& (x .> 0) .& (y .> 0) .& (x .< uplim)
    return x[ps] .+ offset, y[ps]
end

"""
    measure the msd from genetic lineages
"""
function lineageMSD!(msd, counts, gx, gy, srcs, phylo, pbc; offset = 2)
    for src in srcs
        yref = gy[src]
        xref = gx[src]
        current = src; check = true
        @inbounds while current != 0
            if check && abs(gy[current] - gy[src]) <= offset
                xref = gx[current]
                yref = gy[current]
                current = phylo[current]
                continue
            else
                check = false
            end
            x = pbc(abs(xref - gx[current])) ^2
            y = abs(yref-gy[current]) + 1
            msd[y] += x
            counts[y] += 1
            current = phylo[current]
        end
    end
end;

"""
    measure the msd from the shortest paths
"""
function spMSD!(msd, counts, gx, gy, tracks,  pbc; nl=1, offset = 2)
    for (src, trcks) in tracks
        yref = gy[src]
        xref = gx[src]
        for i in eachindex(trcks)
            i <= nl || continue
            trck = trcks[i]
            check = true
            for j in eachindex(trck)
                if check && abs(gy[trck[j]] - gy[src]) <= offset
                    xref = gx[trck[j]]
                    yref = gy[trck[j]]
                    continue
                else
                    check = false
                end
                x = pbc(abs(xref - gx[trck[j]])) ^2
                y = abs(yref-gy[trck[j]]) + 1
                msd[y] += x
                counts[y] += 1
            end
        end
    end 
end

function processSamples(partition; nl = 1, offset = 2)
    params = rand(partition)
    width, height = params.dims
    pvals = []

    # build map geometry
    gxy, _ = buildMap(params.dims, gNodes(), NamedTuple, (n,x,y)->(x = (x - 0.5*(y%2)), y = y))

    gx = getfield.(gxy, :x); gy = getfield.(gxy, :y)

    # periodic boundaries
    pbc = periodicBoundaries(div(width, 2), width)

    # lineage msd
    msd = zeros(Float64, height + 1)
    counts = zeros(Float64, height + 1)

    # SP msd
    spmsd = zeros(Float64, height + 1)
    spcounts = zeros(Float64, height + 1)

    for ithp in partition
        push!(pvals, getfield(ithp, Symbol(ithp.parameter)))

        # measure tracks msd
        if params.density > 0.0 && islink(ithp.path * "/spaths.jld2")
            tracks = load(readlink(ithp.path * "/spaths.jld2"), "tracks")
            spMSD!(spmsd, spcounts, gx, gy, tracks, pbc; nl = nl, offset = offset)
        end

        # measure lineage msd
        jldopen(joinpath(ithp.path, "data.jld2"), "r") do file
            for trial in fContains(keys(file), "trial")
                srcs = file[trial]["sources"]
                phylo = file[trial]["phylo"]

                lineageMSD!(msd, counts, gx, gy, srcs, phylo, pbc; offset = offset)
            end
        end
    end
    return (msd ./ counts), (spmsd ./ spcounts), pvals
end;

function psub(basePath, fSet, offset; overwrite=true)
    rep = first(fSet) # representative configuration

    saveName = joinpath(basePath, join(splitpath(rep.path)[[end - 2, end]], "/")) |> mkpath
    saveName *=  "/msd_processed_data.jld2"

    (isfile(saveName) && ~overwrite) && return saveName

    msd, spmsd, pvals = processSamples(fSet; offset = offset)

    jldsave(saveName; 
        msd = msd,
        spmsd = spmsd,
        config = rep,
        configs = fSet,
        pvalue = measurement(mean_and_std(pvals)...),
    )
    return saveName
end

function picker(path, idxs::Vector{Int} = Int[])
    if isempty(idxs)
        foreach(x->println(x[1],") ",x[2]), enumerate(readdir(path)))
    else
        return readdir(path, join = true)[idxs]
    end
    return nothing
end

# ****************************** main section

# analysis inputs
inputs = (nl = 5, offset = 1, ow = false)

basePath = readline("analysis/lineageMSD/output/defaults/basedirectory.txt") |> mkpath
sampleSet =  (@isdefined pths) ? pths : readlines("analysis/lineageMSD/output/defaults/defaultSample.txt")

fSets = append!(map(getCollections, sampleSet)...)

# process fSets
nameSet = fetch.([@spawn psub($basePath, $fSet, $(inputs.offset); overwrite = inputs.ow) for fSet in fSets])

# (linear) fit model
@. model(x, p) = p[1] * x + p[2]

relabel = Dict("intensity" => L"I", "density"=> L"\phi")

# lineage MSD vs front position
set_theme!(theme_minimal(), fontsize=20) # ColorScheme

let
    fig = Figure(resolution=(700,600), backgroundcolor=:white)

    ax = Axis(fig[1,1], 
        xlabel=L"l_f",
        ylabel="Laterial MSD", 
        xscale=log10, yscale=log10,
        xlabelsize = 38, ylabelsize = 28,
        backgroundcolor=:white,
    )

    # hexagonal lattice scaling
    hexscl = (1 , 1)
    offset = inputs.offset

    xrange = LinRange(5E1, 0.75E3, 10)
    l1 = lines!(ax, xrange, (130 .* xrange) .^ (3/3), color=:black, linewidth = 3)
    text!(ax, L"\alpha=1", position = (140, 2.7E4), align = (:left, :center), rotation = 0.50, fontsize = 28)

    xrange = LinRange(3.5E1, 0.75E3, 10)
    l1 = lines!(ax, xrange, (1 .* xrange) .^ (4/3), color=:black, linewidth = 3)
    text!(ax, L"\alpha=4/3", position = (100, 0.7E3), align = (:left, :center), rotation = 0.60, fontsize = 28)
    
    # xrange = LinRange(10, 1E3, 10)
    # xrange = LinRange(10, 1E3, 10)
    # l2 = lines!(ax, xrange, (1 .* xrange) .^(3/3), color=:black)
    # text!(ax, L"\alpha=1", position = (200, 290), align = (:left, :center), rotation = 0.55, fontsize = 28, linewidth = 5)

    # scanned parameter
    parameter = ""
    # mn, mx = Inf, -Inf
    pvalues = []

    cmap = colorscheme_alpha(roma, 0.7)

    # for i in eachindex(nameSet)
    # foreach(t -> println("$(t[1])  $(t[2])"), enumerate(nameSet))
    for i in [1, 3, 4, 8]
    # for i in [1, 6, 7, 8]
    # for i in [1, 2, 3]
        # load data
        fname = nameSet[i]
        dm = load(fname)

        # simulation configuration
        config = dm["config"]
        parameter = config.parameter
        label = string(round(value(dm["pvalue"]), digits=2))

        col = get(cmap, i / length(fSets)) # colorize

        # remove NaN and Inf
        x, y = hexscl .* passFltr(dm["msd"];
            uplim = 900, offset = 0
        )
        lines!(ax, x, y , color=col, linewidth = 3)

        if config.density > 0.0 && config.intensity > 0
            x, y = hexscl .* passFltr(dm["spmsd"]; 
                uplim = 900, offset = 0
            )
            lines!(ax, x, y, color=col, linewidth = 2.5, linestyle=:dash)
            # lines!(ax, x, y, color=(:red, 0.35), linewidth = 2.5, linestyle=:dash)
        end
        push!(pvalues, getfield(config, Symbol(config.parameter)))
        # display(fig)
        # empty!(ax)
    end

    # axislegend(parameter, position=:lt, orientation=:horizontal, nbanks=4)

    # ranges for pvalues
    mn, mx = extrema(pvalues)

    Colorbar(fig, bbox=BBox(400,650,100,150), colormap=cmap, colorrange=(mn,mx), flipaxis = true, label = relabel[parameter], vertical=false, ticks=0:0.3:mx, labelsize=35)

    ylims!(0.5,10^5)

    display(fig)
    save("analysis/lineageMSD/output/figures/lineagemsd_$(pths |> first |> splitdir |> last).png", fig)
end;