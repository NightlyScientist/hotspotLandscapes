import JLD2: jldsave
import Colors: RGBA
import ColorSchemes: ColorScheme, brg, get, roma, amp

function colorscheme_alpha(
    cscheme::ColorScheme,
    alpha::T = 0.5;
    ncolors = 12,
) where {T<:Real}
    return ColorScheme([
        RGBA(get(cscheme, k), alpha) for k in range(0, 1, length = ncolors)
    ])
end

function survivalFreqImg!(fig, dm)
    ax = Axis(fig[2, 2], backgroundcolor = :transparent)
    rowsize!(fig.layout, 1, Relative(9 / 10))
    heatmap!(ax, 
        reshape(repeat(dm.freq_ancestors,2), (length(dm.freq_ancestors), 2)), lowclip=:white
    )

    Colorbar(fig[3, 2]; label = "Survival Probability", 
        ticksize=20, tickalign=1, ticklabelsize = 22,
        colormap = amp,
        vertical = false, flipaxis = false, labelsize = 35
    )

    hidedecorations!(ax; grid = true)
    xlims!(ax, (1,length(dm.freq_ancestors)))
    rowgap!(fig.layout, 3)
end

function trailsImg(env, lhmap, trails, dm)

    fig = Figure(resolution = (500, 500) .* (√(3), 1.5), backgroundcolor = :transparent)
    ax = Axis(fig[1, 2])

    # axis syling
    hidedecorations!(ax; grid = true)
    limits!(ax, 1, dm.config.dims[1], 1, dm.config.dims[2])

    # construct heatmap from lineages
    h = heatmap!(lhmap, colorrange = (3, 600), colormap=amp, lowclip=:white)

    Colorbar(fig[1, 1], h; label = "Lineage Visitation Frequency",
        width = 10, ticksize=20, ticklabelsize=20, tickalign = 1, labelsize=28, flipaxis = false
    )

    colgap!(fig.layout, 3)

    # show shortest paths
    heatmap!(ax, reshape(trails, dm.config.dims), colorrange = (0.4, 0.5),
        highclip = (:blue, 0.5)
    )

    # the hotspots
    heatmap!(ax, reshape(env, dm.config.dims), colorrange = (1.5, 1.6),
        highclip = (:black, 0.35)
    )

    # label hotspots by visitation status
    # if ~isnothing(preds.pHS) && haskey(preds, :pHS) && ~isempty(preds.pHS.spHS)

    #     cmap = ColorScheme([RGBA(1,1,1, 0.0), RGBA(0.4,0.8,0.0,0.5), RGBA(0.5,0.0,0,0.5), RGBA(1,1,0,.5)])

    #     # threshold = floor(Int, cutoff * maximum(lhmap))
    #     threshold = 0.05

    #     hsColors = zeros(Float64, size(htspts))
    #     # hsColors[preds.pHS.spHS .== 1] .= 1/3
    #     # hsColors[preds.pHS.freqHS .>= threshold] .= 2/3
    #     hsColors[(preds.pHS.freqHS .>= threshold) .& (preds.pHS.spHS .== 1)] .= 1
    #     @info preds.pHS.freqHS

    #     # 0: transparent, 1: sps, 2: lineages, 3: both
    #     x = getfield.(htspts, 1)
    #     y = getfield.(htspts, 2)
    #     scatter!(ax, x, y, color=get(cmap, hsColors), colormap=cmap)
    # end

    return fig
end;

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function accImg(dms::Vector{NamedTuple})

    parameterName = first(dms).config.parameter

    fig = Figure(resolution=(800,750))
    ax = Axis(fig[1, 1], 
        ylabel = "Overlap Fraction", xlabel = L"\phi",
        xlabelsize = 38, ylabelsize = 25,
        yticklabelsize = 25, xticklabelsize = 25
    )

    # hotspots -> agreement with geoemtric paths
    x, xerr = getfield.(dms, :parameter) |> unzip
    y, yerr = getfield.(dms, :acc_hotspots) |> unzip

    # draw accuracy graph for hotspots
    scatter!(ax, x, y, markersize = 35, color=get(roma, 1), 
        label = "Hotspots (M)"
    )
    errorbars!(ax, x, y, xerr, color=get(roma, 1), direction=:x)
    errorbars!(ax, x, y, yerr, color=get(roma, 1))

    # estimates of background for hotspot measure
    bg, bgerr = getfield.(dms, :acc_hotspots_bg) |> unzip

    lines!(ax, x, fill(value(first(bg)), length(x)),
        color=:red, 
    )
    text!(ax, "Random Hotspot Selection", position = (2.5, 0.34), align = (:left, :center),rotation = 0.0, fontsize = 23)

    # ancestors -> agreement with geoemtric paths
    y, yerr = getfield.(dms, :acc_ancestors) |> unzip

    # draw accuracy graph for ancestors
    @info foreach(t->println("$(t[1]) ) $(t[2])"), enumerate(collect(zip(x,y))))
    scatter!(ax, x[2:end], y[2:end], color=get(roma, 0),
        label = "Ancestors (K)", markersize = 35
    )
    errorbars!(ax, x[2:end], y[2:end], xerr[2:end], color=:black, direction=:x)
    errorbars!(ax, x[2:end], y[2:end], yerr[2:end], color=get(roma, 0))

    # plot estimate using a "mean captured counts"
    # y = getfield.(dms, :acc_ancestors_auc)
    # scatter!(ax, x[2:end], y[2:end], color=:green, markersize=25)

    ylims!(ax, (0.2, 1.05))
    axislegend(ax, position=:rb, labelsize=25)
    return fig
end
