include("../base/model.jl")
include("../tools/containers.jl")
include("../tools/parseFiles.jl")

using JLD2, CairoMakie, ColorSchemes, Colors
using ArgParse
import .GraphModel: gNodes, buildMap
import ColorSchemes: tab10, ColorScheme
import Colors: RGBA

function circle(h, k, r, precision=1000)
    θ = LinRange(0, 2π, precision)
    h .+ r .* sin.(θ) .* 1.5, k .+ r .* cos.(θ) .* sqrt(3)
end

function colorscheme_alpha(cscheme::ColorScheme, alpha::T = 0.5; 
        ncolors=12) where T<:Real
    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end


# cmdline = (input="/Users/andromeda/Documents/Projects/Research/lineage-landscapes/storage/simulations/circle/D-0.3271,I-100.0,R-150,XY-(500, 500),HR-(1, 500),trials-1,samples-50", lineages=true, hotspots=true, colors=true)

# cmdline = (input="/Users/andromeda/Documents/Projects/Research/lineage-landscapes/storage/simulations/circle/D-0.003,I-100.0,R-20,XY-(500, 1000),HR-(1, 1000),trials-1,samples-50", lineages=true, hotspots=true, colors=true)

# cmdline = (input="/Users/andromeda/Documents/Projects/Research/lineage-landscapes/storage/simulations/circle/D-0.003,I-7.0,R-20,XY-(500, 1000),HR-(1, 1000),trials-1,samples-50", lineages=true, hotspots=true, colors=true)

# cmdline = (input="/Users/andromeda/Documents/Projects/Research/lineage-landscapes/storage/simulations/circle/D-0.1103,I-9.0,R-87,XY-(500, 500),HR-(1, 500),trials-1,samples-50", lineages=true, hotspots=true, colors=true)
cmdline = (input="/Users/andromeda/Documents/Projects/Research/lineage-landscapes/storage/simulations/animations/D-0.1087,I-4.0,R-150,XY-(500, 500),HR-(1, 500),trials-1,samples-200", lineages=true, hotspots=true, colors=true)

file = jldopen(joinpath(cmdline.input, "data.jld2"), "r")
env = file["env"]
cli = file["cli"]
htspts = file["htspts"]

cli.animation == false && return println("no animation found")

# theme and color
set_theme!(figure_padding = 00)
colors = colorscheme_alpha(ColorSchemes.grays, 0.3, ncolors=3)
greens = colorscheme_alpha(ColorSchemes.Greens, 0.5, ncolors=2)

# create figure
fig = Figure(resolution=cli.dims .* (√(3), 1.5), backgroundcolor=:transparent)
# fig = Figure(resolution=cli.dims, backgroundcolor=:transparent)

ax = Axis(fig[1,1], backgroundcolor=:white)

# set image limits
hidedecorations!(ax)
limits!(ax, 1, cli.dims[1], 1, cli.dims[2])

# used for updating lineages
tmp = zeros(Int64, size(env))

shiftValues(x) = x == 0 ? 0 : 5


# create object for snapshots: list->tuple(growth, lineages)
snapshots = file["trial_1"]["animation"]
close(file)

# i = 41
# i = 49
# i = 43
i = 150
_, ids, lineages = snapshots[i]

# exmpty ax
empty!(ax)

# id colors
if cmdline.colors
    heatmap!(ax, reshape(ids, cli.dims), colormap=:hsv, colorrange=(1,cli.dims[1]), lowclip=:transparent)
else
    heatmap!(ax, shiftValues.(reshape(ids, cli.dims)), colorrange=(1,2), lowclip=:transparent, highclip=(:green, 0.2))
end

# lineages
if cmdline.lineages && length(lineages) > 0
    tmp .= 0; tmp[lineages] .= 3
    color = cmdline.colors ? :black : :black
    heatmap!(ax, reshape(tmp, cli.dims), colorrange=(1,2), lowclip=:transparent, highclip=color)
end

# hotspots
if cmdline.hotspots
    heatmap!(ax, reshape(env, cli.dims), colorrange=(2,3), lowclip=:transparent, colormap=colors)
end

display(fig)


# scale bar
lines!(ax, [10,110], [cli.dims[2] - 10, cli.dims[2] - 10], linewidth=14, color=:black)
text!(ax, L"100", position = (60, cli.dims[2] - 30), align = (:center, :center), fontsize = 40)

display(fig)

# scatter!(ax, [250], [150], markersize=40)
# scatter!(ax, [250], [350], markersize=40)

# add parabla fit
xrange = collect(10:1:490)

# vertex at bottom
# scaleFactor = 1.153333
scaleFactor = 1/1.5
p(x, x_0, I, R, hy) = @. (I + 1) / (4 * I * R) * (x - x_0)^2 + hy - R* I / (I + 1)

yrange = p(xrange, 250, cli.intensity+100, scaleFactor * cli.radius, first(htspts)[2])
lines!(ax, xrange, yrange, linewidth=10, color=(:red, 0.65))

# vertex at center
# p(x, x_0, I, R, hy) = @. (I + 1) / (4 * I * R) * (x - x_0)^2 + hy
# yrange = p(xrange, 250, cli.intensity, cli.radius / 1.5, first(htspts)[2])
# lines!(ax, xrange, yrange, linewidth=10, color=(:black, 0.65))

# text!(ax, L"y = a x^2 + \delta", position = (400, 50), align = (:center, :center), fontsize = 70, color=:white)

display(fig)

# x, y = circle(250,250, cli.radius)
# lines!(ax, x, y, linewidth=10, color=:teal)
# display(fig)

# save("/Users/andromeda/Documents/Projects/Research/lineage-landscapes/storage/images/single_hotspot_snapshot_v3.pdf", fig)