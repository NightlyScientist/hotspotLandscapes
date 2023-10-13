using JLD2
using FileIO
using CairoMakie
using StatsBase
using ColorSchemes
using EasyFit
using ArgParse
import Base.Threads: @spawn, fetch
import Measurements: value, uncertainty, measurement, ±
import EasyFit.LsqFit: curve_fit, stderror
import ColorSchemes: roma, amp, get
import Colors: RGBA

include("../../src/tools/parseFiles.jl");
include("../../src/tools/binUtils.jl");
include("../../src/tools/containers.jl")

function picker(path, idxs::Vector{Int} = Int[])
    if isempty(idxs)
        foreach(x->println(x[1],") ",x[2]), enumerate(readdir(path)))
    else
        return readdir(path, join = true)[idxs]
    end
    return nothing
end

function processSamples(partition; nbins = 20, smpl = true)

    pvals = []
    finalTimes = Float64[]

    jldopen(rand(partition).path * "/data.jld2", "r") do file
        for trial in fContains(keys(file), "trial")
            push!(finalTimes, file[trial]["finalTime"])
        end
    end

    ft = mean(finalTimes)
    f = binParams([0.0, 1.2 * ft]; n = nbins)

    # define binning variables
    fnames = [:sectorsize, :fpos, :fvar, :nroots]
    cnts = []
    binned = []

    for i in eachindex(fnames)
        push!(cnts, zeros(Int64, nbins))
        push!(binned, zeros(Float64, nbins))
    end

    for ithp in partition

        push!(pvals, getfield(ithp, Symbol(ithp.parameter)))
        file = jldopen(joinpath(ithp.path, "data.jld2"), "r")

        for trial in fContains(keys(file), "trial")
            sctr = file[trial]["sectors"]
            frnt = file[trial]["front"]

            # time
            t = getfield.(sctr, 1)

            # mean sector size
            smoothing!(binned[1], cnts[1], t, getfield.(sctr, 2), f)

            # number of lineage roots (unique population ids)
            smoothing!(binned[4], cnts[4], t, getfield.(sctr, 4), f)

            # mean front position
            smoothing!(binned[2], cnts[2], t, getfield.(frnt, 2), f)

            # variance about mean front position
            smoothing!(binned[3], cnts[3], t, getfield.(frnt, 3), f)
        end
        close(file)
    end

    time = (ft / nbins) .* (collect(1:nbins) .- 0.5)

    measures = Dict{Symbol,Any}()
    for i in eachindex(fnames)
        data = (binned[i] ./ cnts[i])
        measures[fnames[i]] = (time = time, data = data)
    end
    return measures, pvals
end;

function colorscheme_alpha(
    cscheme::ColorScheme,
    alpha::T = 0.5;
    ncolors = 12,
) where {T<:Real}
    return ColorScheme([
        RGBA(get(cscheme, k), alpha) for k in range(0, 1, length = ncolors)
    ])
end

# parallel run and save data
function prun(ptn, basedir, nbins)
    fname = splitdir(first(ptn).path) |> last
    measures, pvals = processSamples(ptn; nbins = nbins)

    jldsave(joinpath(basedir, fname * ".jld2"); 
        mean_and_std = mean_and_std,
        measures = measures,
        pvals = pvals,
        config = ptn |> first,
        configs = ptn
    )
    return pvals
end

# convert (phi, r) -> λ := n^-2
lambda(ϕ, r) = sqrt(-π * r^2 / log(1 - ϕ))

# base project director
basedir = "/storage/jgonzaleznunez/lineage-landscapes/scratch/analysis/frontMotion" |> mkpath

# find data fSets from path
inputs =
    join(readlines(joinpath(pwd(), "analysis/front", "path.input")), " ") |>
    Meta.parse |>
    eval |>
    namedtuple
    
pth = inputs.path
fSets = getCollections(pth)

# process fSets
pvals = fetch.([@spawn prun($ptn, $basedir, $(inputs.nbins)) for ptn in fSets])

# ColorScheme
set_theme!(theme_minimal(), fontsize = 28)
update_theme!(
    Theme(
        Axis = (
            xgridvisible = false,
            ygridvisible = false,
            xtrimspine = false,
            ytrimspine = (false, true),
        ),
    ),
)

insetPos(x, y, width, aspect) = (x, x + width, y, y + width * aspect);

relabel = Dict("intensity" => L"I", "density"=> L"\phi")

@. model(x, p) = p[1] * x + p[2]

# plots
plotOne = true

let prm = string(fSets[1][1].parameter)
    fig = Figure(resolution = (800, 800), backgroundcolor=:transparent, figure_padding = 0)

    if ~plotOne
        inset_1 = Axis(
            fig,
            bbox = BBox(insetPos(541, 190, 250, 0.8)...),
            # yscale = log10,
            xlabel = relabel[prm],
            ylabel = L"\Lambda",
            xgridvisible = false,
            ygridvisible = false,
            xtrimspine = false,
            ytrimspine = (false, true),
            # xticks = 0:0.25:0.8,
            backgroundcolor=:transparent,
            xlabelsize = 37,
            ylabelsize = 40,
        )
    end

    inset_2 = Axis(
        fig,
        bbox = BBox(insetPos(210, 450, 250, 1)...),
        ylabel = "Front Speed",
        xlabel = relabel[prm],
        xgridvisible = false,
        ygridvisible = false,
        xtrimspine = true,
        ytrimspine = (false, true),
        xticks = WilkinsonTicks(4),
        xlabelsize = 40,
        backgroundcolor=:transparent,
    )

    ax = Axis(
        fig[1, 1],
        yscale = log10,
        xscale = log10,
        xlabel = L"\bar{h}",
        ylabel = L"l_s(\bar{h})",
        xgridvisible = false,
        ygridvisible = false,
        xlabelsize = 45,
        ylabelsize = 50,
        backgroundcolor=:transparent
    )

    scatterplot = []

    cmap = colorscheme_alpha(roma, 0.9; ncolors = length(fSets))
    cmap_2 = colorscheme_alpha(roma, 0.6; ncolors = length(fSets))

    for i in eachindex(fSets[sortperm(pvals, by=mean)])
        # i == lastindex(pvals) && continue

        fname = splitdir(first(fSets[i]).path) |> last
        dm = load(joinpath(basedir, fname * ".jld2"))

        config = dm["config"]
        measures = dm["measures"]

        # color points by scanned parameter
        pval = measurement(mean_and_std(pvals[i])...)
        label = string(round(value(pval), digits = 2))
        col = get(cmap, i / length(fSets))

        # hotspot separation length scale
        λ = config.density > 0. ? lambda(config.density, config.radius) : Inf

        # draw length scale vs parameter in inset
        plotOne || scatter!(inset_1, value(pval), λ + 2*config.radius, color = col, markersize=19)

        # extract abscissa and ordinate for number of roots
        x = measures[:nroots].time
        y = (config.dims[1] ./ measures[:nroots].data)

        # extract front position
        ht = measures[:fpos].data

        # estimate front propogation speed from linear fit
        fit = fitlinear(x[ht.>10], ht[ht.>10])
        standardErr = curve_fit(model, x[ht.>10], ht[ht.>10], [0.5,0.5]) |> stderror

        v₀ = 1

        # draw front speed and errors
        errorbars!(inset_2, [value(pval)], [fit.a ./ v₀], [uncertainty(pval)]; direction = :x,)
        errorbars!(inset_2, [value(pval)], [fit.a ./ v₀], [standardErr[1]], direction = :y)
        scatter!(inset_2, value(pval), fit.a ./ v₀, color = col, markersize=25)

        # draw average sector size (main data)
        ps = (ht .!== NaN) .& (y .!== NaN) .& (ht .> 0) .& (y .> 0) .& (x .< 1000)
        lines!(ax, ht[ps][1:end-10], y[ps][1:end-10], color = col, linewidth = 5)

        if i == 1
            xrange = collect(1E2:1.2E3)
            lines!(ax, xrange, 3.1.* xrange .^ (2/3), color=:black)

            sfit = fitlinear(log.(ht[ps][10:end-10]), log.(y[ps][10:end-10]))
            @info sfit

            text!(ax, "α = 2/3", position = (200, 140), align = (:left, :center),rotation = 0.72, fontsize = 28)
        end

        # add mark to indicate location of λ + 2R
        frst = findfirst(y[ps] .>= λ + 20)
        if frst !== nothing
            push!(scatterplot, (ht[ps][frst], y[ps][frst], get(cmap_2, i / length(fSets))))
            # push!(scatterplot, (ht[ps][frst], y[ps][frst], col))
        end
    end

    mn, mx = extrema(reduce(vcat, pvals))
    Colorbar(fig, bbox=BBox(250,520,110,130), colormap=cmap,  flipaxis = true, vertical = false, colorrange=(mn,mx), label = relabel[prm], labelsize=40, ticks = 0:0.5:mx)

    prm == "density" && scatter!(ax, 
        getfield.(scatterplot, 1), 
        getfield.(scatterplot, 2), 
        color = getfield.(scatterplot, 3),
        marker='|',  markersize=55
    )

    display(fig)
    fname = splitpath(pth) |> last
    @info fname
    save("analysis/front/output/figures/sector_sizes_$(fname).png", fig)
end;