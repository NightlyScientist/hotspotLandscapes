import ColorSchemes: get, ColorScheme, hsv
using CairoMakie, Colors
using DataStructures

function colorscheme_alpha(cscheme::ColorScheme, alpha::T = 0.5; 
        ncolors=12) where T<:Real
    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

function landVisual(dims, env, colors, phylogeny = nothing; drawScatter=false)
    # create figure and axis
    fig = Figure(resolution = (500, 500), backgroundcolor = :transparent)
    ax = Axis(fig[1, 1], backgroundcolor = :transparent)

    # visual representation of growth
    drawSnapshot!(
        ax,
        dims,
        env;
        lineages = collect(keys(phylogeny)),
        colors = colors,
        drawScatter = drawScatter
    )

    return fig, ax
end

function convrt(index, w) 
    x = mod(index, 1:w)
    y = fld(index - x, w) + 1
    return x, y
end

# plot sectors, lineages, hotspots
function drawSnapshot!(ax, dims, env; lineages = nothing, colors = nothing, drawScatter=false)

    alphaColors = colorscheme_alpha(hsv, 1; ncolors=dims[1])

    # colors (ids)
    if !isnothing(colors)
        heatmap!(
            ax,
            reshape(colors, dims),
            # colormap = :darktest,
            colormap = alphaColors,
            colorrange = (1, dims[1]),
            lowclip = :transparent,
        )
    end

    # set image limits
    hidedecorations!(ax)
    limits!(ax, 1, dims[1], 1, dims[2])

    # hotspots
    heatmap!(
        ax,
        reshape(env, dims),
        colorrange = (1.5, 1.6),
        lowclip = :transparent,
        highclip = (:black, 0.45),
    )

    # plot lineages as spatial traces
    if ~isnothing(lineages) 
        if drawScatter
            x = Int[]
            y = Int[]
            for lidx in lineages
                tmp_x, tmp_y = convrt(lidx, dims[1])
                push!(x, tmp_x)
                push!(y, tmp_y)
            end
            scatter!(ax, x, y, color=:lightgray, markersize=3)
        else
            tmp = zeros(Int64, size(env))
            tmp[lineages] .= 3
            heatmap!(
                ax,
                reshape(tmp, dims),
                colorrange = (1, 2),
                lowclip = :transparent,
                highclip = :white,
                # highclip = :black,
            )
        end
    end
end
