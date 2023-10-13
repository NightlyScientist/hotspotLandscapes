import JLD2: jldopen, jldsave
import FileIO: load
import StatsBase: mean, kurtosis, std, median, mode, skewness, mean_and_std
import Colors: RGBA
import ColorSchemes: get, amp, roma, hot, ColorScheme, grays, nipy_spectral
import Base.Threads: @spawn, fetch
import Measurements: value, uncertainty, ±, measurement
import DataStructures: DefaultDict
using CairoMakie

include("../../src/tools/parseFiles.jl");
include("../../src/tools/containers.jl")
include("../../src/tools/binUtils.jl")
include("measures.jl")

function colorscheme_alpha(
    cscheme::ColorScheme,
    alpha::T = 0.5;
    ncolors = 12,
) where {T<:Real}
    return ColorScheme([
        RGBA(get(cscheme, k), alpha) for k in range(0, 1, length = ncolors)
    ])
end

unzip(a::Vector{NamedTuple}) = map(x -> getfield.(a, x), fieldnames(eltype(a)))
unzip(a::Vector{Tuple{Any}}) = map(x -> getfield.(a, x), eachindex(a))

function picker(path, idxs::Vector{Int} = Int[])
    if isempty(idxs)
        foreach(x->println(x[1],") ",x[2]), enumerate(readdir(path)))
    else
        return readdir(path, join = true)[idxs]
    end
    return nothing
end

function traverse!(target::Vector{<:Real}, srcs, phylo)
    for src in srcs
        current = src
        @inbounds while current != 0
            target[current] += 1
            current = phylo[current]
        end
    end
end

function lineageHM(data, dims)
    lhmap = zeros(Int64, reduce(*, dims))

    for trial in fContains(keys(data), "trial")
        srcs = data[trial]["sources"]
        phylo = data[trial]["phylo"]
        traverse!(lhmap, srcs, phylo)
    end
    return lhmap
end

function survivingAncestors(fSet)
    survivalnumbers = Int64[]
    for rep in fSet
        jldopen(joinpath(rep.path, "data.jld2"), "r") do data
            for trial in fContains(keys(data), "trial")
                pids = length(unique(data[trial]["pIDs"]))
                push!(survivalnumbers, pids)
            end
        end
    end
    return measurement(mean_and_std(survivalnumbers)...)
end

# parallel run and save data
function psub(basePath, fSet; overwrite = false)
    rep = first(fSet) # representative, collecition-1
    saveName = joinpath(basePath, join(splitpath(rep.path)[[end - 2, end]], "/")) |> mkpath

    jldopen(joinpath(rep.path, "data.jld2"), "r") do input
        if ~isfile("$saveName/lhmap_data.jld2") || overwrite
            lhmap = lineageHM(input, rep.dims)
            jldsave("$saveName/lhmap_data.jld2"; config = rep, lhmap = lhmap, fSet = fSet)
        end
        if ~isfile("$saveName/sa_data.jld2") || overwrite
            sa = survivingAncestors(fSet)
            jldsave("$saveName/sa_data.jld2"; config = rep, sa = sa, fSet = fSet)
        end
    end
end

# **************** main

# analysis inputs
inputs = (nl = 15,  showimg = true, smpl = true, ow = false)

basePath = mkpath(readline("analysis/shortestpaths/output/defaults/basedirectory.txt"))

# collect simulation ensemble and run in parallel
sampleSet =  (@isdefined pths) ? pths : readlines("analysis/shortestpaths/output/defaults/defaultsample.txt")
fSets = append!(map(getCollections, sampleSet)...)
wait.([@spawn psub($basePath, $fSet; overwrite = inputs.ow) for fSet in fSets])

# ********** plotting

# * sample lineage Visitation freq heatmap, chosen at slctr index
let opt_1 = 2, opt_2 = 4
    set_theme!(theme_minimal(), fontsize = 28, figure_padding = 10)
    # fig = Figure(resolution = (800, 800), backgroundcolor = :transparent)
    fig = Figure(resolution = (600, 600) .* (√(3), 1.5), backgroundcolor = :transparent)

    ax = Axis(fig[1, 1], backgroundcolor = :white)
    colors = colorscheme_alpha(grays, 0.3, ncolors = 3)

    name_1 = picker(basePath, [opt_1]) |> first
    name_2 = picker(name_1, [opt_2]) |> first

    @info name_1 name_2

    db = load(joinpath(basePath, name_2, "lhmap_data.jld2"))

    repConfig = db["config"]
    lhmap = reshape(db["lhmap"], repConfig.dims)

    nrmlztn = repConfig.numberTrials * repConfig.dims[1]

    mn = minimum(filter( >(0.0), lhmap ./ nrmlztn))
    colorrange = (mn, 0.8*maximum(lhmap ./ nrmlztn))

    heatmap!(ax, lhmap ./ nrmlztn,
        colormap = amp,
        colorrange = colorrange,
        lowclip = :white,
    )

    Colorbar(fig[2, 1];
        # label = "Lineage Visitation Frequency",
        label = "Survival Probability",
        ticksize = 15,
        tickalign = 1,
        colormap = amp,
        # height = 770,
        vertical = false,
        flipaxis = false,
        colorrange = colorrange,
        labelsize = 40
    )
    rowgap!(fig.layout, 2)

    if repConfig.intensity > 0
        env = load(joinpath(repConfig.path, "data.jld2"), "env")
        heatmap!(ax,
            reshape(env, repConfig.dims),
            colorrange = (2, 3),
            lowclip = :transparent,
            colormap = colors,
        )
    end

    hidedecorations!(ax; grid = true)
    colgap!(fig.layout, 2)
    rowsize!(fig.layout, 1, Relative(1.0))
    display(fig)

    @info splitdir(name_2) |> last
    # save("analysis/shortestpaths/output/figures/siteVisitationHM_$(last(splitdir(name_2))).png", fig)
    save("analysis/shortestpaths/output/figures/siteVisitationHM_sample_2.png", fig)
end

# * histogram of Visitation frequencies at some height
let slctr = 25
    set_theme!(theme_minimal(), fontsize = 28)
    fig = Figure(resolution = (800, 800), backgroundcolor = :white)
    ax = Axis(
        fig[1, 1],
        backgroundcolor = :white,
        xlabel = "Lineage Visitation Frequency",
        ylabel = "Probability",
    )

    # open data file if it exists
    fname = nameSet[slctr] * "/lhmap_data.jld2"
    lhmap, repConfig = load(joinpath(basePath, fname), "lhmap", "config")
    lhmap = reshape(lhmap, repConfig.dims)

    # get lineage heat map
    if haskey(load(dname), "fSet")
        fSet = load(dname, "fSet")
        vals = getfield.(fSet, Symbol(rep.parameter))
        parameter = measurement(mean_and_std(vals)...)
    else
        parameter = measurement(getfield(rep, Symbol(rep.parameter)), 0)
    end

    # filter out zeroes
    x = filter(>(0.0), pinning(reshape(lhmap, repConfig.dims), 500, 10))

    # bin according to x, ignore y, and retrieve bin counts
    x_binned, _, bin_counts = binning(x, copy(x), LinRange(0, maximum(x), 70) |> collect)

    # remove zero-size bins (to use in log scale)
    fltr = (x_binned .> 0) .& (bin_counts .> 0)
    scatter!(ax, x_binned[fltr], bin_counts[fltr] ./ sum(bin_counts))

    display(fig)
    # save(
    #     "analysis/shortestpaths/output/vf_$(ptn.parameter)_$(Measurements.value(parameter)).png",
    #     fig,
    # )
end

# * create animation of survival frequencies, retrieve stats
let idxs = [12]
    set_theme!(theme_light(), fontsize = 28, figure_padding = 1)
    fig = Figure(resolution = (800, 700), backgroundcolor = :transparent)
    ax = Axis(fig[1, 1], backgroundcolor = :white)
    hidedecorations!(ax)

    # hotspot color
    colors = colorscheme_alpha(grays, 0.3, ncolors = 3)

    top_dir = picker(basePath, idxs) |> first |> splitdir |> last

    record(
        fig,
        "analysis/shortestpaths/output/videos/pinningeffect_$(top_dir).mp4";
        framerate = 2,
    ) do io

        # iterate through paths
        for dbSet in picker(basePath, idxs)
            pvals = getfield.(map(x->load(x * "/lhmap_data.jld2", "config"), readdir(dbSet, join=true)), :intensity)

            for fname in readdir(dbSet, join = true)[sortperm(pvals)]

                dname = fname * "/lhmap_data.jld2"
                isfile(dname) || continue

                # open data file if it exists
                lhmap, rep = load(dname, "lhmap", "config")

                # empty previous content from axis
                empty!(ax)

                y₀ = fld(rep.dims[2], 2)
                x = filter(>(0.0), pinning(reshape(lhmap, rep.dims), y₀, rep.radius))

                heatmap!(
                    ax,
                    reshape(lhmap, rep.dims),
                    colormap = amp,
                    colorrange = (1, maximum(lhmap)),
                    lowclip = :transparent,
                )

                env = load(joinpath(rep.path, "data.jld2"), "env")
                heatmap!(
                    ax,
                    reshape(env, rep.dims),
                    colorrange = (2, 3),
                    lowclip = :transparent,
                    colormap = colors,
                )
                recordframe!(io)
            end
        end
    end
end

# * draw stats (e.g.. kurtosis, mean, etc)
let idxs = [13]
    set_theme!(theme_minimal(), fontsize = 28)
    colors = colorscheme_alpha(nipy_spectral, 0.8, ncolors = 6)
    fig = Figure(backgroundcolor = :white, resolution=(800,600))

    ax = Axis(
        fig[1, 1],
        xlabel = "Intensity",
        ylabel = "Lineage Site-Visitation Frequency",
        backgroundcolor = :white,
        rightspinevisible = false,
        xlabelsize = 30
    )

    twin = Axis(
        fig[1, 1],
        xlabel = "Intensity",
        backgroundcolor = :white,
        yaxisposition = :right,
        yticklabelcolor = colors[5],
        rightspinecolor = colors[5],
        ytickcolor = colors[5],
        xlabelsize = 30
    )

    twin_2 = Axis(
        fig[1, 1],
        xlabel = "Intensity",
        backgroundcolor = :white,
        yaxisposition = :right,
        yticklabelcolor = colors[6],
        rightspinecolor = colors[6],
        ytickcolor = colors[6],
        xlabelsize = 30
    )

    # normalize to the base case
    nrm = 1 ± 0

    statsSet = Dict{Tuple,Vector{NamedTuple}}()
    phaseParams = [:radius, :density, :intensity]

    # iterate through paths
    for dbSet in picker(basePath, idxs), fname in readdir(dbSet, join = true)
        dname = joinpath(basePath, fname * "/lhmap_data.jld2")
        isfile(dname) || continue

        lhmap, rep = load(dname, "lhmap", "config")

        if haskey(load(dname), "fSet")
            fSet = load(dname, "fSet")
            vals = getfield.(fSet, Symbol(rep.parameter))
            parameter = measurement(mean_and_std(vals)...)
        else
            parameter = measurement(getfield(rep, Symbol(rep.parameter)), 0)
        end

        lhmap = reshape(lhmap, rep.dims)

        x = filter(>(0.0), pinning(lhmap, fld(rep.dims[2], 2), rep.radius))

        p = map(x->getfield(rep, x), filter(!=(Symbol(rep.parameter)), phaseParams))
        k = (p..., rep.dims...)

        if ~haskey(statsSet, k)
            statsSet[k] = NamedTuple[]
        end

        push!(statsSet[k],
            (
                parameter = parameter,
                mn = mean(x),
                sd = std(x),
                mdn = median(x),
                md = mode(x),
                skns = skewness(x),
                kurt = kurtosis(x),
                rep.intensity,
                rep.radius,
                rep.density
            ),
        )
    end

    for (k, stats) in statsSet
        x = value.(getfield.(stats, :parameter))
        lines!(ax, x, getfield.(stats, :mn), color = colors[1], label = "mean",
            linestyle = :dash, linewidth = 5)
        lines!(ax, x, getfield.(stats, :mdn), color = colors[2], label = "median",
            linestyle = :dash, linewidth = 5)
        lines!(ax, x, getfield.(stats, :md), color = colors[3], label = "mode",
            linestyle = :dash, linewidth = 5)
        lines!(
            ax,
            x,
            getfield.(stats, :sd),
            color = colors[4],
            label = "Standard Deviation",
            linestyle = :dash,
            linewidth = 5,
        )
        lines!(
            twin,
            x,
            getfield.(stats, :skns),
            color = colors[5],
            label = "skewness",
            linestyle = :dash,
            linewidth = 5,
        )
        lines!(
            twin_2,
            x,
            getfield.(stats, :kurt),
            color = colors[6],
            label = "kurtosis",
            linestyle = :dash,
            linewidth = 5,
        )
        @info k
        break
    end

    Legend(fig, ax, bbox = BBox(330, 350, 270, 300))
    Legend(fig, twin, bbox = BBox(550, 570, 410, 450))
    Legend(fig, twin_2, bbox = BBox(550, 570, 290, 310))

    ylims!(ax, 0.0, 0.2)
    ylims!(twin_2, 0.0, 13)
    # xlims!(ax, 0.0, 4)
    # xlims!(twin, 0.0, 4)
    # xlims!(twin_2, 0.0, 4)

    display(fig)

    top_dir = pths |> first |> splitdir |> last
    # save("analysis/shortestpaths/output/stats_$(top_dir).pdf", fig)
end

# * overlay surviving_ancestors from multiple fSets
let idxs = collect(11:19)
# let idxs = [12, 13, 14, 15, 11]
# let idxs = [5,6,7,8,9,2]
    # pinPrm(ϕ, I) = @. √(ϕ)*(I / (I + 1))
    pinPrm(ϕ, I) = @. √(-log(1 - ϕ))*(I / (I + 1))

    set_theme!(theme_minimal(), fontsize = 28)
    fig = Figure(backgroundcolor = :transparent, resolution=(800,800))

    ax = Axis(
        fig[1, 1],
        backgroundcolor = :transparent,
        ylabel = "Normalized Number of Surviving Ancestors",
        xlabel = L"\phi^ {1/2} \frac{I}{I+1}",
        # xlabel = L"I",
    )

    phaseParams = [:radius, :density, :intensity]
    statsSet = Dict{String,Vector{NamedTuple}}()

    for dbSet in picker(basePath, idxs), fname in readdir(dbSet, join = true)
        dname = fname * "/sa_data.jld2"
        isfile(dname) || continue

        sa, rep = load(dname, "sa", "config")
        if haskey(load(dname), "fSet")
            fSet = load(dname, "fSet")
            vals = getfield.(fSet, Symbol(rep.parameter))
            parameter = measurement(mean_and_std(vals)...)
        else
            parameter = measurement(getfield(rep, Symbol(rep.parameter)), 0)
        end

        k = splitdir(dbSet) |> last

        ~haskey(statsSet, k) && (statsSet[k] = NamedTuple[])

        push!(statsSet[k],
            (
                parameter = parameter,
                val = value(sa),
                err = uncertainty(sa),
                lx = rep.dims[1],
                ly = rep.dims[2],
                density = rep.density,
                radius = rep.radius,
                intensity = rep.intensity,
                prmName = rep.parameter
            )
        )
    end

    parameterName = :density
    # mn, mx = value.(getfield.(values(statsSet), parameterName)) |> extrema
    paramValues = []
    mn, mx = Inf, 0.41

    cmap = colorscheme_alpha(roma, 0.6)

    @info "number in statsSet: " length(statsSet)
    for (k, stats) in statsSet
        x, y, yerr, lx, ly, ϕ, R, I, pn =  map(x -> getfield.(stats, x), keys(first(stats)))

        color = get(cmap, mean(ϕ) / mx)
        if mn >= mean(ϕ)
            mn = mean(ϕ)
        end
        # color = :black

        # x = @. (I / (I + 1))
        x = @. sqrt(ϕ) * (I / (I + 1))
        # x = @. sqrt(sqrt(ϕ)) * (I / (I + 1))
        # x = @. sqrt(-0.5*log(1-ϕ)) * (I / (I + 1))
        # x = @. ϕ^(1/4) * (I / (I + 1))
        # x = @. sqrt(-log(1 - ϕ)) * (I / (I + 1))
        x = x[sortperm(I)]
        y = y[sortperm(I)]

        scatter!(ax, value.(x), y ./ maximum(y) , color = color, markersize=35)
        lines!(ax, value.(x), y  ./ maximum(y), color = color)
        errorbars!(ax, value.(x), y  ./ maximum(y), yerr ./ maximum(y) , color = color)

        # lines!(ax, collect(0:0.1:8), pinPrm.(mean(ϕ), collect(0:0.1:8)), color = color)
        # find when sqrt(phi) I / (i+1) reaches 50% of its value
        # hlines!(ax, pinPrm(mean(ϕ), 1), color = color)

        @info extrema(ϕ), extrema(R), extrema(I)
        # push!(paramValues, mean(ϕ))
    end 

    # text!(ax, L"$\sqrt{\phi} \frac{I}{I+1}$", position = (2, 0.1), align = (:left, :center), rotation = 0.05, fontsize = 20, color = :red)

    # @info mn, mx
    # Colorbar(fig, bbox=BBox(100,460,130,150), colormap=roma,  flipaxis = true, vertical = false, colorrange=(mn,mx), label = L"ϕ", labelsize=45)


    ylims!(ax, 0.1,1.1)
    # xlims!(ax, -.01, .3)
    display(fig)
    # save("analysis/shortestpaths/output/figures/surviving_ancestors.png", fig)
    save("/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/analysis/shortestpaths/output/figures/scaled_survAnc.png", fig)
end