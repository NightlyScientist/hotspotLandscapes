include("../base/environment.jl")
include("../base/model.jl")

using CairoMakie
using ColorSchemes
using Colors
import .Model: gNodes, buildMap, predefinedModels, hexGraph, simulate!, linearStart!, resetGraph!

function colorscheme_alpha(cscheme::ColorScheme, alpha::T = 0.5; 
        ncolors=12) where T<:Real
    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

# global parameters
dimensions = 20, 20
regionBounds = 1, dimensions[2]
density = 0.1 
intensity = 4
radius = 1
stopline = 17

# image path for figures
imagePath = joinpath(pwd(), "storage/images/hexLattice")
mkpath(imagePath)

# construct data model containers
dataModels = predefinedModels()

# build hex map
xy, cntns = hexGraph(dimensions)

# reset models
Model.resetModels!(dataModels)

# reset graph properties
resetGraph!(xy, cntns)

env, _ = uniformHotspots(dimensions, density, radius, regionBounds)

# xy, cntns = buildMap(
#     dimensions,
#     gNodes(),
#     Tuple,
#     (n, x, y) -> (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1)),
# )

# gx = getfield.(xy, 1)
# gy = getfield.(xy, 2)

gx = getfield.(xy, :x) ./ (sqrt(3))
gy = ((getfield.(xy, :y) .- 1) ./ 1.5 ) .+ 1

# perform mini simulation
opts = (radius=radius, density=density, intensity=intensity, dims=dimensions, numberSamples=50,)

# initial population
active = linearStart!(xy, cntns, dimensions, env)

# do simulation
final = simulate!(xy, cntns, active, env, opts, dataModels; stopline=stopline)

sources = final.fSnapshot

# determine the phylo of the periphery population 
phylo, branches = Model.genealogy(xy, sources)

# set_theme!(theme_light())
set_theme!(theme_minimal())

cmap = let
    precmap = to_color.(Makie.PlotUtils.palette(:darktest))

    popat!(precmap, 6)
    precmap

    tmp = pop!(precmap)

    pushfirst!(precmap, tmp)

    cmap = ColorScheme(precmap)
end

# cylic colormap
cmap = ColorSchemes.hsv
alphaColors = colorscheme_alpha(ColorSchemes.hsv, 1; ncolors=opts.dims[1])

# let canvasLX = 450, ms = cld(canvasLX, 15)
let canvasLX = 450, ms = canvasLX / 16
    fig = Figure(resolution = (canvasLX, canvasLX * sqrt(3) / 2),
        backgroundcolor = :transparent
    )
    ax = Axis(fig[1, 1],
        # xticks = 1:1:10
    )

    # grid as hexagons (black outline)
    scatter!(ax, gx, gy, marker = :hexagon, markersize = ms + 0.5, color = :black)

    # grid as hexagons (white fill)
    scatter!(ax, gx, gy, marker = :hexagon, markersize = ms - 0.5, color = (:white, 1.0))

    # color initial population
    color = collect(1:dimensions[1])
    scatter!(
        ax,
        gx[color],
        gy[color],
        marker = :hexagon,
        markersize = ms - 10,
        color=color[color],
        colormap = alphaColors,
        # colormap = :binary,
    )

    #  color population
    fltr = getfield.(xy, :filled) .== true
    fltr[color] .= false
    color = getfield.(xy, :ID_2)
    scatter!(
        ax,
        gx[fltr],
        gy[fltr],
        marker = :hexagon,
        markersize = ms - 0.5,
        color=color[fltr],
        colormap = alphaColors,
        # colormap = :binary,
    )

    # hotspots
    scatter!(
        ax,
        gx[env.==2],
        gy[env.==2],
        marker = :hexagon,
        markersize = ms - 1,
        color = (:black, 0.4),
    )

    # lineages
    # lineages = keys(phylo) |> collect
    # scatter!(
    #     ax,
    #     gx[lineages],
    #     gy[lineages],
    #     marker = :hexagon,
    #     markersize = ms/4,
    #     color = (:black, 1),
    # )

    for source in sources
        current = source
        idxs = Int64[]

        while current != 0
            push!(idxs, current)
            current = phylo[current]
        end

        lines!(
            ax,
            gx[idxs],
            gy[idxs],
            color = (:white, 1),
            linewidth = 5
        )
    end

    # lineages but only the front
    scatter!(
        ax,
        gx[final.fSnapshot],
        gy[final.fSnapshot],
        marker = :hexagon,
        markersize = ms/3,
        color = (:black, 1),
    )

    # hidedecorations!(ax)
    save(imagePath * "/hexlattice.png", fig)
    display(fig)
end