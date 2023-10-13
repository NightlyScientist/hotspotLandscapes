import Base.Threads: @spawn, fetch
import JLD2: jldopen, jldsave
import FileIO: load
import StatsBase: mean, kurtosis, std, median, mode
import Colors: RGBA
import ColorSchemes: get, amp, roma, hot, ColorScheme, brg
import Measurements: value, uncertainty, ±, measurement
import EasyFit.LsqFit: curve_fit, stderror
import EasyFit: fitlinear
using CairoMakie
using DataStructures
using LaTeXStrings

include("../../src/base/model.jl")
include("../../src/tools/binUtils.jl")
include("../../src/tools/parseFiles.jl")
include("../../src/tools/trees.jl")
include("../../src/tools/containers.jl")

import .GraphModel: gNodes, buildMap

function picker(path, idxs::Vector{Int} = Int[])
    if isempty(idxs)
        foreach(x->println(x[1],") ",x[2]), enumerate(readdir(path)))
    else
        return readdir(path, join = true)[idxs]
    end
    return nothing
end

function colorscheme_alpha(
    cscheme::ColorScheme,
    alpha::T = 0.5;
    ncolors = 12,
) where {T<:Real}
    return ColorScheme([
        RGBA(get(cscheme, k), alpha) for k in range(0, 1, length = ncolors)
    ])
end

# convert (phi, r) -> λ := n^-2
lambda(ϕ, r) = sqrt(π * r^2 / (-log(1 - ϕ)));

# binning done with: https://arxiv.org/pdf/1207.5578.pdf
function processSamples(partition; cutoff = 0.1, nbins = 20, smpl = false)
    # get random finalTime
    rp = rand(partition)
    pvals = []

    ft = let finalTimes = Float64[]
        # estimate time binning from from sample of trials
        jldopen(joinpath(rp.path, "data.jld2"), "r") do file
            for trial in fContains(keys(file), "trial")
                push!(finalTimes, file[trial]["finalTime"])
            end
        end
        mean(finalTimes)
    end

    # graph positions
    xy, _ = buildMap(
        rp.dims,
        gNodes(),
        NamedTuple,
        # dev: x axis should be shifted, but doesn't matter -> y = const
        (n, x, y) -> (x = x, y = 1 + 1.5 * (y - 1)),
    )

    # define binning variables
    width = rp.dims[1]
    mdist = ceil(Int64, min(cutoff * width, width / 2))

    jtx = Vector{Vector{Float64}}(undef, 1 + mdist)
    foreach(i -> jtx[i] = zeros(Int64, nbins), eachindex(jtx))

    # avoid fintie size effects via P.B.C.
    pbc = periodicBoundaries(floor(Int64, width / 2), width)
    f = binParams([0, ft]; n = nbins)

    gx = getfield.(xy, :x)

    for ithp in partition
        push!(pvals, getfield(ithp, Symbol(ithp.parameter)))

        file = jldopen(joinpath(ithp.path, "data.jld2"), "r")

        for (i, trial) in enumerate(fContains(keys(file), "trial"))
            tree = file[trial]["tree"]
            pids = file[trial]["pIDs"]
            srcs = file[trial]["sources"]

            tjtx = Jτx(mdist, pbc, tree, srcs, pids, gx)

            # bin all data into J(τ, Δx₀)
            foreach(i -> updateHist!(jtx[i], tjtx[i], f), eachindex(tjtx))
            smpl && i > 50 && break
        end
        close(file)
        smpl && break
    end
    return (jtx=jtx, mft=ft, τ=(ft / nbins) .* (collect(1:nbins) .- 0.5)), pvals
end

# using P.B.C get the arc distance between two points
function periodicBoundaries(halfwidth::Int, width::Int)
    f = let h = halfwidth, w = width
        x -> (x > h) ? w - x : x
    end
    return f
end

function Jτx(mdist, pbc, tree, srcs, pids, gx)
    jtx = Vector{Vector{Float64}}(undef, 1 + mdist)
    foreach(i -> jtx[i] = Vector{Float64}(), eachindex(jtx))

    a::Int = 0
    b::Int = 0
    for i in eachindex(srcs), j = i+1:length(srcs)
        a = srcs[i]
        b = srcs[j]

        dx::Int = pbc(abs(gx[a] - gx[b]))
        0 < dx <= mdist || continue

        p = lca_times(tree, a, b, pids[i], pids[j])
        p == 0 && continue

        τ = (tree[a].time + tree[b].time - 2 * tree[p].time) * 0.5

        # skip 
        iszero(τ) && continue

        # fill Jτx
        push!(jtx[dx], τ)
    end
    return jtx
end

function tauHist(jtx, times, maxdist)
    xdata = Float64[]

    for i in eachindex(jtx)
        Float64(i) > maxdist && continue

        # skip if there are no values for this separation
        sum(jtx[i]) > 0 || continue

        # remove entries with zero counts
        flt = jtx[i] .> 0
        xdata = vcat(xdata, (jtx[i].*times)[flt])
    end

    p = sortperm(vcat(xdata, times))
    bins = findall(>(length(xdata)), p)
    bin_count = Int.(diff(bins) .- 1)

    return bin_count, times
end

function T2_xt(jtx, times, maxdist)
    t2x = NTuple{2,Float64}[]
    # @info length(jtx)
    for i in eachindex(jtx)
        i <= maxdist || continue
        sum(jtx[i]) > 1 || continue
        push!(t2x, (i, sum((jtx[i] .* times)) / sum(jtx[i])))
    end
    return getfield.(t2x,1), getfield.(t2x,2)
end;

function rescaleJtx(jtx, times, nbins, maxdist)
    xdata = Float64[]
    ydata = Float64[]

    for i in eachindex(jtx)
        Float64(i) > maxdist && continue

        # skip if there are no values for this separation
        sum(jtx[i]) > 1 || continue

        # we scale by sqrt(3) to get ditsances on hex grid
        # scl = (sqrt(3) * Float64(i)) ^ (3/2)
        scl = Float64(i)^(3 / 2)

        # normalize by the integral
        # nrml = sum(jtx[i] .* (first(diff(times))))
        nrml = sum(jtx[i])

        # remove entries with zero counts
        flt = jtx[i] .> 0

        # append data to containers
        xdata = vcat(xdata, times[flt] ./ scl)
        ydata = vcat(ydata, jtx[i][flt] .* scl ./ nrml)
    end

    # binEdges = LinRange(extrema(xdata)..., nbins) |> collect
    binEdges = logspace(log.(2, extrema(xdata))..., nbins; base = 2.0) |> collect
    x_binned, y_binned, bin_counts = binning(xdata, ydata, binEdges)

    flt = bin_counts .> 0
    return x_binned[flt], y_binned[flt]
end

function psub(basePath, fSet, nbins, cutoff; smpl = false, overwrite = true)
    rep = first(fSet) # representative, collecition-1
    saveName = joinpath(basePath, join(splitpath(rep.path)[[end - 2, end]], "/")) |> mkpath

    # fname = splitdir(first(fSet).path) |> last
    # fname = basedir * "/$(fname)_nbins-$(nbins)_cutoff_$(cutoff).jld2"

    (isfile(saveName * "/coalescences.jld2") && ~overwrite) && return

    d, pvals = processSamples(fSet; nbins = nbins, cutoff = cutoff, smpl = smpl)

    jldsave(saveName * "/coalescences.jld2"; 
        jtx = d.jtx,
        tau = d.τ,
        mft = d.mft,
        config = fSet |> first,
        configs = fSet,
        nbins = nbins,
        pval = measurement(mean_and_std(pvals)...),
        pvals = pvals
    )
end

unzip(a::Vector{NamedTuple}) = map(x->getfield.(a, x), fieldnames(eltype(a)))
unzip(a::Vector{Tuple{Any}}) = map(x->getfield.(a, x), eachindex(a))

# ****************************** main section

# analysis inputs
inputs = (cutoff = 0.2, ow = true, smpl = false, nbins = 4000)

basePath = readline("analysis/commonAncestry/output/defaults/basedirectory.txt") |> mkpath
sampleSet =  (@isdefined pths) ? pths : readlines("analysis/commonAncestry/output/defaults/defaultsample.txt")

fSets = append!(map(getCollections, sampleSet)...)

wait.([@spawn psub($basePath, $fSet, $(inputs.nbins), $(inputs.cutoff); overwrite = inputs.ow, smpl = inputs.smpl) for fSet in fSets])

# where to place insets
insetPos(x, y, width, aspect) = (x, x + width, y, y + width * aspect);

# (linear) fit model
@. model(x, p) = p[1] * x + p[2]


# coalescence rate plot
let nbins = 1400, maxdists = (250,10), idxs = [3,4]
    set_theme!(theme_minimal(), fontsize = 35)
    fig = Figure(; resolution = (1000, 600), backgroundcolor=:transparent)
    cmap = roma

    # colsize!(fig.layout, 1, Relative(2/3))
    ax11 = Axis(
        fig[1, 1],
        xlabel = L"\tau / \Delta x_0^{3/2}",
        ylabel = L"\Delta x_0^{3/2} \cdot J(\tau \; | \; \Delta x_0)",
        xscale = log10,
        yscale = log10,
        xlabelsize = 45,
        ylabelsize = 45
    )

    fig_inset = Figure(backgroundcolor=:transparent)
    ax22 = Axis(fig_inset[1,1], 
        # xlabel = string(inputs.parameter), 
        # ylabel = L"γ",
        backgroundcolor=:transparent,
        xlabelsize = 50,
        ylabelsize = 70,
        xticks = [0, 0.6]
    )

    ax22twin = Axis(fig_inset[1,1], 
        yaxisposition = :right,
        yticklabelcolor = :black,
        rightspinecolor = :black,
        ytickcolor = :black,
        # ylabel="β",
        backgroundcolor=:transparent,
        xlabelsize = 50,
        ylabelsize = 50,
        xticks = [0, 0.6]
    )

    # extra paths to include
    appendPths = picker(basePath, idxs)
    fnames = append!(readdir.([appendPths...], join=true)...) .* "/coalescences.jld2"

    mxpvals = value.(map(t->load(t, "pval"), fnames))
    betas = []

    # for dbSet in picker(basePath, idxs), fname in readdir(dbSet, join = true)
    for i in eachindex(sortperm(mxpvals))

        # load data
        dm = load(joinpath(fnames[i]))

        # simulation configuration
        config = dm["config"]

        parameter = config.parameter
        pval = dm["pval"]

        label = string(round(value(pval), digits=2))

        # mean and std of values grouped
        mval = value(pval)
        vstd = uncertainty(pval)

        # jtx and associated variables
        jtx = dm["jtx"]
        τ = dm["tau"]

        # color
        col = get(cmap, value(pval) / maximum(mxpvals))

        # rescale coalescence rate
        rt, rj = rescaleJtx(jtx, τ, nbins, maxdists[1])

        # normalize axes
        nrj = (rj./ 1)[rt.<100.0]
        trt = rt[rt.<100.0]

        # find power law for τ/x >> 1 
        fltr = (3.1 * 10^0 .< trt .< 9.0 * 10^1)
        fit = fitlinear(log.(trt[fltr]), log.(nrj[fltr]))

        # find linear fit
        fit2 = curve_fit(model, log.(trt[fltr]), log.(nrj[fltr]), [0.5,0.5])

        # confidence intervals
        confIntvl = stderror(fit2)

        # plot alpha
        scatter!(ax22, [mval], [fit.a], color=(:gray, 0.8), markersize=36)

        # plot uncertainty of alpha
        errorbars!(ax22, [mval], [fit.a], [confIntvl[1]];
            whiskerwidth = 12,
            color = :black,
            direction = :y,
        )
        errorbars!(ax22, [mval], [fit.a], [uncertainty(pval)];
            whiskerwidth = 12,
            color = :black,
            direction = :x,
        )

        # rescale factor, c, appearing in (c τ) / Δx
        # rsclng = 1
        rsclng = exp((fit.b ± confIntvl[2]) / (fit.a ± confIntvl[1])) 

        # draw rescaled coalescence for τ / x >> 1, with rscling shift
        scatter!(ax11, value(rsclng) .* trt, nrj, markersize = 10, color = col)

        rscl = 1 ./ rsclng

        # scatter!(ax22twin, [mval], [fit.b], color=:dodgerblue)

        scatter!(ax22twin, [value(pval)], [value(rscl)], color=:green, marker=:dtriangle, markersize=35)

        # add beta to list
        push!(betas, value(rscl))

        errorbars!(ax22twin, [value(pval)], [value(rscl)], 
            [uncertainty(rscl)];
            whiskerwidth = 12,
            direction = :y,
            color = :black
        )

        # lines!(
        #     ax11,
        #     trt[fltr],
        #     exp(fit.b) .* trt[fltr] .^ fit.a,
        #     markersize = 5,
        #     color = :blue,
        #     label = "$(round(fit.a,digits=2))",
        # )
    end

    # dev: calculate λ(ϕ)
    mn, mx = filter(!iszero, mxpvals) |> extrema
    x = LinRange(mn, mx, 100)
    y = lambda.(x, 10) .+ 10 
    y = y .* ( maximum(betas) / maximum(y))
    lines!(ax22twin, x .- 0.02, y, color=:blue, linewidth=5, label=L"λ(\phi, R)")

    # dev: add color bar to fit_inset
    # Colorbar(fig_inset[1,1], label = L"\phi", vertical = false, height=15, ticksize=10, tickalign=1, flipaxis=false, colorrange=(0., maximum(mxpvals)), colormap=roma, bbox=BBox(001,100,100,100))

    Colorbar(fig_inset[2,1], label = L"\phi", vertical = false,
    width = Relative(4 / 4), height=15, ticksize=10, tickalign=1, flipaxis=false, colorrange=(0., maximum(mxpvals)), colormap=roma, labelsize=45)

    hidespines!(ax22, :t, :r, :b)
    hidespines!(ax22twin, :t, :r, :b)

    # limits!(ax11, (5 * 10^-2, 6 * 10^2), (1 * 10^-5, 8 * 10^0))
    ylims!(ax22, (-1.3, -1.9))
    axislegend(ax22twin, position=:rt, labelsize=45, linewidth=80, markersize=80)
    display(fig)
    # Legend()
    # display(fig_inset)
    save("analysis/commonAncestry/output/ca_scaled.png", fig)
    # save("analysis/commonAncestry/output/ca_unscaled.png", fig)
    save("analysis/commonAncestry/output/ca_inset.png", fig_inset)
end

# * T2 measures
let
    appendPths = picker(basePath, [1,2,3,4,5])

    cmap = roma
    set_theme!(theme_minimal(), fontsize = 28)
    fig = Figure(; resolution = (830, 600), backgroundcolor=:white)

    # colsize!(fig.layout, 1, Relative(2/3))
    ax = Axis(fig[1, 1],
        xlabel = L"\phi",
        ylabel = L"T_2/T_2^{\mathrm{uniform}}",
        # ylabel = "Mean Common Ancestry Time",
        xlabelsize = 40,
        ylabelsize = 35
    )

    N₀ = 1

    epths = readdir("/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/storage/scratch/analysis/CommAncs", join=true)
    
    xvals = Float64[]
    yvals = Float64[]

    parameterValues = []
    for dbSet in epths
        mxpvals = maximum(value.(map(t->load(t, "pval"), readdir(dbSet, join=true))))
        append!(parameterValues, value.(mxpvals))
    end

    mn, mx = extrema(parameterValues)

    append!(epths, appendPths)

    dataSets = Dict()

    for dbSet in epths
    # for dbSet in picker(basePath, idxs)
        # mxpvals = maximum(value.(map(t->load(t * "/coalescences.jld2", "pval"), readdir(dbSet, join=true))))
        # mxpvals = maximum(value.(map(t->load(t, "pval"), readdir(dbSet, join=true))))

        for fname in readdir(dbSet, join = true)
            if contains(fname, ".jld2")
                dm = load(joinpath(fname))
            else
                dm = load(joinpath(fname * "/coalescences.jld2"))
            end
            repConfig = dm["configs"] |> first
            jtx = dm["jtx"]
            τ = dm["tau"]
            mft = dm["mft"]
            pval = dm["pval"]

            if repConfig.intensity > 3 && value(pval) < 0.2
                @info "pvalue " pval
            end

            # color = get(roma, value(pval) / mx)

            distances, t2xt = T2_xt(jtx, τ, 55)

            value(pval) == 0. && (N₀ = mean(t2xt) / mft)

            y = measurement(mean_and_std(t2xt)...) / mft / N₀

            if value(pval) > 0.
                # scatter!(ax, 
                #     value(pval), value(y),
                #     markersize = 30,
                #     marker=:utriangle,
                #     color = color
                # )

                k = repConfig.intensity

                if haskey(dataSets, k)
                    push!(dataSets[k], (value(pval), value(y)))
                else
                    dataSets[k] = [(value(pval), value(y))]
                end
            end

            vlines!(ax, 1 - 0.68)
            text!(ax, L"1 - \phi'", position = (0.3, 0.61), align = (:left, :center), rotation = π/2, fontsize = 35)

            # errorbars!(ax, 
            #     value(pval), value(y), 
            #     [uncertainty(y)];
            #     whiskerwidth = 12,
            #     direction = :y,
            #     color = :black
            # )
        end
    end

    cols = colorscheme_alpha(brg, 0.4; ncolors = 10)

    for col in sort(keys(dataSets) |> collect)
        ds = dataSets[col]
        x = getfield.(ds, 1)
        y = getfield.(ds, 2)
        x_bin, y_bin, b_c = binning(x, y, 0.0:0.00005:0.8)

        @info x_bin[b_c .> 0]
        scatter!(ax, x_bin[b_c .> 0], y_bin[b_c .> 0], 
            markersize = 30,
            # marker=:circle,
            color = get(cols, col / maximum(keys(dataSets))),
            # color = cols[col],
            label = "I = $(col)"
        )
    end

    # srtr = sortperm(xvals)
    # xvals = xvals[srtr]
    # yvals = yvals[srtr]
    # fit = fitlinear(log10.(xvals[0 .< xvals .< 0.3]), log10.(yvals[0 .< xvals .< 0.3]))
    # xr =  collect(0.02:0.001:0.3)
    # lines!(ax, xr, exp(fit.b).*(xr .^(fit.a)), color=:black )
    # Colorbar(fig, bbox=BBox(220,620,650,670), colormap=roma, colorrange=(mn,mx), flipaxis = false, label = "ϕ", labelsize = 35, vertical=false)
    # ylims!(ax, 0.5, 0.8)
    axislegend(ax, position = :rb)
    display(fig)
    save("analysis/commonAncestry/output/figures/mean_t2_v4.png", fig)
end