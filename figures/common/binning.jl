"""
    bin(x::Vector{<:Number}, y::Array{<:Number}, xbin_edges::Union{Vector{<:Number},StepRangeLen{}}; method::F = mean) where {F<:Function}

Fast binning of `y` as a function of `x`. Returns `x_binned`, `y_binned`, `bin_count`. `method` 
specifies approach used for aggregating binned values. NOTE: `x_binned` == 0 and 
`y_binned` == 0 when `bin_count` == 0.

# see https://discourse.julialang.org/t/performance-optimization-of-a-custom-binning-funciton/91616
"""
function binning(
  x::Vector{<:Number}, y::Array{<:Number}, xbin_edges::Union{Vector{<:Number},StepRangeLen{}}; method::F=mean
) where {F<:Function}

  # find bin breakpoints
  p = sortperm(vcat(x, xbin_edges))
  bins = findall(>(length(x)), p)

  # initialize outputs
  bin_count = Int.(diff(bins) .- 1)
  x_binned = zeros(eltype(x), length(bin_count))
  y_binned = zeros(eltype(y), length(bin_count), size(y, 2))

  # calculate binned metrics
  @inbounds for i in findall(bin_count .> 0)
    x_binned[i] = method(@view(x[p[(bins[i] + 1):(bins[i + 1] - 1)]]))
    y_binned[i] = method(@view(y[p[(bins[i] + 1):(bins[i + 1] - 1)]]))
  end
  return x_binned, y_binned, bin_count
end

logspace(start, finish, num::Int64; base::Float64=10.0) = base .^ LinRange(start, finish, num)

function rebin(x, y, nbins, method=mean)
  xbin_edges = collect(LinRange(extrema(x)..., nbins))
  return binning(x, y, xbin_edges; method=method)
end

# https://stats.stackexchange.com/questions/74843/binning-by-equal-width
function binMapping(x; n=2, a=1)
  (mn, mx) = extrema(x)
  r = mx - mn
  h = r / (n - a)
  δ = (n * h - r) / 2
  f = let mn = mn, δ = δ, h = h
    t -> ceil(Int64, (t - mn + δ) / h)
  end
  return f
end

function updateHist!(counts, x, f::Function)
  mx = length(counts)
  @inbounds for n in eachindex(x)
    k = f(x[n])
    if 0 < k <= mx
      counts[k] += 1
    end
  end
end

function smoothing!(binned, counts, x, y, f::Function)
  mx = length(counts)
  @inbounds for n in eachindex(x)
    k = f(x[n])
    # investigate: should we clamp() instead?
    if 0 < k <= mx
      counts[k] += 1
      binned[k] += y[n]
    end
  end
end

function groupSame(x, y)
  N = [[t] for t in y]
  w = fill(true, length(x))
  for i in 2:axes(x)
    if x[i - 1] == x[i] && w[i - 1]
      append!(N[i], N[i - 1])
      w[i - 1] = false
    end
  end
  return N[w]
end