include("../../src/common/indexTools.jl")
import Measurements: value, uncertainty, ±, measurement, Measurement
import StatsBase: mean, kurtosis, std, median, mode, skewness, mean_and_std

module FilterTools
mask(haystack::Vector{String}, needle::String) = filter(t -> contains(t, needle), haystack)
mask(haystack::Base.KeySet, needle::String) = filter(t -> contains(t, needle), collect(haystack))

function mask(x::Vector{<:Real}, y::Vector{<:Real}; xlbound=-Inf, xhbound=Inf, scale=identity)
  _x = scale.(x)
  _y = scale.(y)
  _mask = isfinite.(_x) .& isfinite.(_y) .& (_x .> scale(xlbound)) .& (_x .< scale(xhbound))
  return x[_mask], y[_mask]
end
end

module FitModels
import LsqFit: curve_fit, stderror, coef
export LinearFitModel, linearfit

struct LinearFitModel
  a::Float64
  b::Float64
  a_err::Float64
  b_err::Float64
  sumSqRes::Float64
end

@. linearModel(x, p) = p[1] * x + p[2]

function linearfit(x, y; xlbound=-Inf, xhbound=Inf, scale=identity)::Union{Missing,LinearFitModel}
  _x, _y = mask(x, y; xlbound=xlbound, xhbound=xhbound, scale=scale)
  isempty(_x) && isempty(_y) && return missing
  fit = curve_fit(linearModel, scale.(_x), scale.(_y), [0.0, 1.5])
  errors = sqrt.(stderror(fit))
  sumsqres = sum(abs2, scale.(_y) .- linearModel(scale.(_x), coef(fit)))
  return FitModels.LinearFitModel(coef(fit)..., errors..., sumsqres)
end

function mask(x::Vector{<:Real}, y::Vector{<:Real}; xlbound=-Inf, xhbound=Inf, scale=identity)
  _x = scale.(x)
  _y = scale.(y)
  _mask = isfinite.(_x) .& isfinite.(_y) .& (_x .> scale(xlbound)) .& (_x .< scale(xhbound))
  return x[_mask], y[_mask]
end
end

namedtuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

pinningParameter(ϕ, I) = @. √(-log(1 - ϕ)) * (I / (I + 1))
pinningParameter(ϕ, I, scalefactor::Function=sqrt) = @. (I / (I + 1)) * scalefactor(ϕ)

lambda_sqr(ϕ, r) = sqrt(π * r^2 / (-log(1 - ϕ)))
lambda_hex(ϕ, r) = sqrt(-π * (2 / sqrt(3) * r)^2 / log(1 - ϕ))
lambda(ϕ, r, func::Function) = ϕ > 0.0 ? func(ϕ, r) : Inf

periodicBoundaries(w::Int, x::Int) = (x > cld(w, 2)) ? w - x : x
periodicBoundaries(w::Int, x::Float64) = (x > /(w, 2)) ? w - x : x

module DataFrameTools
using DataStructures, CSV, DataFrames, DataFramesMeta, StatsBase
export Results, convert, unzip, bounds, groupframe

mutable struct Results{T}
  x::Union{Missing,Vector{T}}
  y::Union{Missing,Vector{Any}}
  colors::Union{Missing,Vector{Any}}
end

_replace(symbol::Symbol, old::Symbol=:density, new::Symbol=:group) = first(replace([symbol], old => new))
convert(df::DataFrame, opt=:group) = groupby(df, _replace(opt))
convert(df::SubDataFrame, opt=:group) = groupby(df, replace(opt))
convert(df::GroupedDataFrame, opt=:group) = df

unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))
unzip(a::Vector{NamedTuple}) = map(x -> getfield.(a, x), fieldnames(eltype(a)))
unzip(a::Vector{Tuple{Any}}) = map(x -> getfield.(a, x), eachindex(a))

bounds(df::DataFrame, opt) = extrema(df[!, opt])
bounds(df::SubDataFrame, opt) = extrema(df[!, opt])

function bounds(df::GroupedDataFrame, opt)
  _values = []
  foreach(subframe -> push!(_values, extrema(subframe[!, opt])...), df)
  return extrema(_values)
end

function groupframe(df, values::Vector{<:Real}, index::Int; opts=(:density, :intensity))
  first_opt = _replace(first(opts), :density, :group)
  subframe_1 = @subset(df, $(first_opt) .== values[index])
  # subframe_1 = @subset(df, $(first(opts)) .== values[index])

  subdivide = mean(subframe_1[!, opts[1]]) * mean(subframe_1[!, opts[2]]) > 0.0
  subdivide || return groupby(subframe_1, _replace(last(opts), :density, :group))
  return subframe_1
end
end
