include("../../src/base/model.jl")
include("../../src/base/containers.jl")
include("../../src/base/environment.jl")

using JLD2, CairoMakie, ColorSchemes, Colors
using ArgParse
import .Model: gNodes, buildMap
import ColorSchemes: tab10, ColorScheme
import Colors: RGBA

function circle(h, k, r, precision=1000)
  θ = LinRange(0, 2π, precision)
  h .+ r .* sin.(θ) .* 1.5, k .+ r .* cos.(θ) .* sqrt(3)
end

function colorscheme_alpha(cscheme::ColorScheme, alpha::T=0.5; ncolors=12) where {T<:Real}
  return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1; length=ncolors)])
end

storage_path = "/storage/disorderedLandscapes/simulations/animated"
cmdline = (
  input = "",
  lineages=true,
  hotspots=true,
  colors=true
)

density, radius, width, height = load(joinpath(cmdline.input, "Opts.jld2"), "density", "radius", "width", "height")
intensity = load(joinpath(cmdline.input, "Opts.jld2"), "intensity")
animation = load(joinpath(cmdline.input, "Opts.jld2"), "animate") 
objs = load(joinpath(cmdline.input, "objects.jld2"), "objs")
env, _ =  applyObstacles!(objs, radius, width, height)

# env = file["env"]
# cli = file["cli"]
# htspts = file["htspts"]
htspts = objs

animation == false && return println("no animation found")

# theme and color
set_theme!(; figure_padding=00)
colors = colorscheme_alpha(ColorSchemes.grays, 0.3; ncolors=3)
greens = colorscheme_alpha(ColorSchemes.Greens, 0.5; ncolors=2)

# create figure
fig = Figure(; resolution=(width, height) .* (√(3), 1.5), backgroundcolor=:transparent)
# fig = Figure(resolution=cli.dims, backgroundcolor=:transparent)

ax = Axis(fig[1, 1]; backgroundcolor=:white)

# set image limits
hidedecorations!(ax)
limits!(ax, 1, width, 1, height)

# used for updating lineages
tmp = zeros(Int64, size(env))

shiftValues(x) = x == 0 ? 0 : 5

# create object for snapshots: list->tuple(growth, lineages)
extras = load(joinpath(cmdline.input, "data_extras.jld2"))
snapshots = extras["trial_1/animation"]

i = 85
_, ids, lineages = snapshots[i]

let
# exmpty ax
empty!(ax)

# id colors
if cmdline.colors
  heatmap!(ax, reshape(ids, (width, height)); colormap=:hsv, colorrange=(1, width), lowclip=:transparent)
else
  heatmap!(ax, shiftValues.(reshape(ids, cli.dims)); colorrange=(1, 2), lowclip=:transparent, highclip=(:green, 0.2))
end

# lineages
if cmdline.lineages && length(lineages) > 0
  tmp .= 0
  tmp[lineages] .= 3
  color = cmdline.colors ? :black : :black
  heatmap!(ax, reshape(tmp, (width, height)); colorrange=(1, 2), lowclip=:transparent, highclip=color)
end

# hotspots
if cmdline.hotspots
  heatmap!(ax, reshape(env, (width, height)); colorrange=(2, 3), lowclip=:transparent, colormap=colors)
end

display(fig)

# scale bar
cli = (dims = (width, height), intensity = intensity, radius = radius,)
lines!(ax, [10, 110], [cli.dims[2] - 10, cli.dims[2] - 10]; linewidth=14, color=:black)
text!(ax, L"100"; position=(60, cli.dims[2] - 30), align=(:center, :center), fontsize=40)

display(fig)

xrange = collect(10:1:490)
scaleFactor = 1 / 1.153333
p(x, x_0, I, R, hy, r) = @. (I + 1) / (4 * I * R) * (x - x_0)^2 + hy - 1.15*r * I / (I + 1)

yrange = p(xrange, 250, cli.intensity,  scaleFactor * cli.radius, first(htspts)[2], radius)
lines!(ax, xrange, yrange; linewidth=10, color=(:gray, 0.95))

display(fig)
end

_path = "Projects/disorderedLandscapes/workspace/images/single_circle_snapshot/fig_1B.pdf"
save(_path, fig)
