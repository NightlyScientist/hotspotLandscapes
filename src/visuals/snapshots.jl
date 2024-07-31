import ColorSchemes: get, ColorScheme, hsv
include("../common/indexTools.jl")
using CairoMakie, Colors
using DataStructures

function alphaColor(cscheme::ColorScheme, alpha::T=0.5; ncolors=12) where {T<:Real}
  return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1; length=ncolors)])
end

function landscape(graph, dims, data_path=nothing)
  width, height = dims
  fig, _ = heatmap(
    reshape(getfield.(graph, :ID_1), (width, height)); colorrange=(1.1, 1.2), lowclip=:white, highclip=:black
  )
  save(joinpath(data_path, "ID_1.png"), fig)
end

function landVisual(dims, graph, phylogeny=nothing; data_path=nothing, drawScatter=false)
  fig = Figure(; size=(500, 0.87 * 500 * dims[2] / dims[1]), backgroundcolor=:transparent)
  ax = Axis(fig[1, 1]; backgroundcolor=:transparent)

  env = getfield.(graph, :ID_1)
  colors = getfield.(graph, :ID_2)

  drawSnapshot!(ax, dims, env; lineages=phylogeny, colors=colors, drawScatter=drawScatter)

  isnothing(data_path) || save(joinpath(data_path, "ID_2.png"), fig)
  return fig, ax
end

unzip(tuple) = map(collect, zip((Tuple(e) for e in tuple)...))

# doc: plot sectors, lineages, hotspots
function drawSnapshot!(ax, dims, env; lineages=nothing, colors=nothing, drawScatter=false)
  hidedecorations!(ax)
  limits!(ax, 1, dims[1], 1, dims[2])

  # .draw genetic labels (ids)
  kwargs = (colormap=alphaColor(hsv, 1; ncolors=dims[1]), colorrange=(1, dims[1]), lowclip=:transparent)
  isnothing(colors) || heatmap!(ax, reshape(colors, dims); kwargs...)

  # .draw hotspots
  heatmap!(ax, reshape(env, dims); colorrange=(1.5, 1.6), lowclip=:transparent, highclip=(:black, 0.45))

  # .draw lineages as spatial traces
  isnothing(lineages) || return nothing
  if drawScatter
    x, y = unzip(reconstructCoordinates(collect(keys(lineages)), dims...))
    scatter!(ax, x, y; color=:lightgray, markersize=3)
  else
    tmp = zeros(Int64, size(env))
    tmp[lineages] .= 3
    heatmap!(ax, reshape(tmp, dims); colorrange=(1, 2), lowclip=:transparent, highclip=:white)
  end
end
