using Revise
using StaticArrays
using StructArrays
using CairoMakie

include("../base/environment.jl")

# two ideas: 
# (1) iterate through all hotspots, each time filling all other hotspot
# as being visited by this one, then next in height setting the ancestor(s)
# as coming from shorstest hotspot within the list of visited ones


mutable struct Hotspot
    x::Float64
    y::Float64
    touched::Bool
    angles::Vector{Float64}
    previous::Vector{Int64}
    Hotspot(x, y) = new(x,y, false, [], [])
end

function _inParabola(x_0, y_0, x, y, θ, a)
    # apply rotaton
    dx = x - x_0
    dy = y - y_0

    # check if in the correct side of the plane
    s = sin(θ); c = cos(θ)

    δy = -dx * s + dy * c
    δy >= 0 || return false

    δx = dx * c + dy * s

    return a * δx ^2 <= δy
end

function _distances(hotspots::Vector{Hotspot})
    N = length(hotspots)
    distances = zeros(Float64, (N,N))

    x = getfield.(hotspots, :x)
    y = getfield.(hotspots, :y)

    for i in axes(hotspots,1)
        @. distances[i,:] = sqrt((x - hotspots[i].x) ^2 + (y - hotspots[i].y) ^2)
    end
    return distances
end

function _findNext(hotspots::Vector{Hotspot}, distances::Matrix, idx, a, θ=0)
    x_0 = hotspots[idx].x
    y_0 = hotspots[idx].y

    inside = (Inf, 0, idx)
    outside = (Inf, 0, idx)

    θ₀ = θ

    for i in eachindex(hotspots)
        idx != i || continue

        θ₀ = atan(hotspots[i].y - y_0, hotspots[i].x - x_0)

        if _inParabola(x_0, y_0, hotspots[i].x, hotspots[i].y, θ, a)
            if distances[idx, i] <= inside[1]
                inside = (distances[idx, i], θ₀, i)
            end
        elseif distances[idx, i] <= outside[1]
            outside = (distances[idx, i], θ₀, i)
        end
    end
    return inside, outside
end

function parabola(a, k, x, width, precision=100)
    xrange = LinRange(0, width, precision) |> collect
    xrange, a .* (xrange .- x).^2 .+ k
end

position(hotspots, loc) = [hotspots[loc[3]].x], [hotspots[loc[3]].y]

cli = (landscape=nothing, hsRegion=nothing, dims=(100,100), radius=2, envType="", density=0.1)

env, htspts, cli = createEnvironment!(cli)

_hotspots = sort(htspts, by=last)
hotspots = [Hotspot(x,y) for (x,y) in zip(first.(_hotspots), last.(_hotspots))]

# pairwise distances
distances = _distances(hotspots)

fig, ax = scatter(getfield.(hotspots, :x), getfield.(hotspots, :y), color=:black)

# check which hotspots are within the parabola
inside, outside = _findNext(hotspots, distances, 1, 0.1)

scatter!(ax, position(hotspots, inside)..., color=:red)
scatter!(ax, position(hotspots, outside)..., color=:blue)

x_0 = hotspots[1].x
y_0 = hotspots[1].y
θ = 0

for i in eachindex(hotspots)
    _inParabola(x_0, y_0, hotspots[i].x, hotspots[i].y, θ, 0.1) |> println
end

xr, yr = parabola(0.1, hotspots[1].y, hotspots[1].x, cli.dims[2])
lines!(ax, xr, yr, linewidth=2)

inside, outside = _findNext(hotspots, distances, inside[3], 0.1, inside[2])
scatter!(ax, position(hotspots, inside)..., color=:red)
scatter!(ax, position(hotspots, outside)..., color=:blue)

# xr, yr = parabola(0.1, hotspots[1].y, hotspots[1].x, cli.dims[2])
# lines!(ax, xr, yr, linewidth=2)

ylims!(ax, 0, cli.dims[2])
xlims!(ax, 0, cli.dims[1])
display(fig)