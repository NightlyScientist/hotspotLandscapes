using JLD2, FileIO, ArgParse, Base.Threads, CSV, DataFrames

module Ancestry
using JLD2, FileIO, Base.Threads
include("../../src/calculations/lineages.jl")
include("../../src/base/trees.jl")
include("../common/theme.jl")
include("../common/binning.jl")
include("../common/modeling.jl")

module AncestryTools
  mutable struct Results
    jtx::Matrix{Float64}
    binSize::Float64
    nbins::Int64
    bins::Vector{Float64}
  end
end

# typemap = Dict("Main.GraphModel.TreeNode" => TreeNode)

# doc: binning done with: https://arxiv.org/pdf/1207.5578.pdf
function mainProcess(dataPath, nbins, cutoffFraction)::AncestryTools.Results
  opts = namedtuple(load(joinpath(dataPath, "Opts.jld2")))
  sources = GenealogicalTree.graphSources(opts)

  # .avoid fintie size effects via P.B.C.
  max_dist = ceil(Int64, min(cutoffFraction * opts.width, opts.width / 2))
  coalRatetx = zeros(Float64, (1 + max_dist, nbins))
  maxTime = 1.0
  binIndex = binMapping([0, maxTime]; n=nbins)

  file = jldopen(joinpath(dataPath, "data_phylo.jld2"), "r")
  trials = LineageTracing.mask(keys(file), "trial")
  for (i, trial) in enumerate(trials)
    geneticLabels::Vector{UInt32} = file[trial]["source_labels"]
    genealogy = file[trial]["phylogeny"]
    branchPoints = file[trial]["branch_points"]
    branchTimes = file[trial]["branch_times"]
    gtree = GenealogicalTree.buildTree(genealogy, geneticLabels, branchPoints, branchTimes, sources)

    if firstindex(trials) == i
      maxTime = file[trial]["end_time"]
      binIndex = binMapping([0, maxTime]; n=nbins)
    end
    Jτx!(coalRatetx, opts.width, max_dist, gtree, sources, geneticLabels, binIndex)
  end
  close(file)

  bin_size = maxTime / nbins
  bins = collect(1:nbins) .* bin_size .- 0.5 * bin_size
  return AncestryTools.Results(coalRatetx, bin_size, nbins, bins)
end

# doc: add coalescence rates to input matrix
function Jτx!(coalRatetx::Matrix{Float64}, width::Int, max_dist::Int, tree, sources, labels::Vector{UInt32}, binIndex)
  a::Int32 = 0
  b::Int32 = 0
  for i in eachindex(sources), j in (i + 1):length(sources)
    a = @inbounds sources[i]
    b = @inbounds sources[j]

    # .assumes that sources are on the same line, a - b == x(a) - x(b)
    dx::Int64 = GenealogicalTree.periodicBoundaries(width, abs(a - b))
    (dx > max_dist || dx == 0) && continue

    p = GenealogicalTree.lca_times(tree, a, b, labels[i], labels[j])
    p == 0 && continue

    τ = @inbounds (tree[a].time + tree[b].time - 2 * tree[p].time) * 0.5

    iszero(τ) && continue
    bin_index::UInt32 = binIndex(τ)
    coalRatetx[dx, bin_index] += 1
  end
end

function rescaleCoalescenceRate(results, maxdist; scale=1.5, nbins=4000)
  _x = Float64[]
  _y = Float64[]

  for dx in axes(results.jtx, 1)
    (Float64(dx) > maxdist || sum(results.jtx[dx, :]) == 0) && continue

    # .scale by sqrt(3) to get distances on hex grid
    # _scale = (sqrt(3) * Float64(dx))^(scale)
    _scale = Float64(dx)^(scale)

    _mask = results.jtx[dx, :] .> 0
    _x = vcat(_x, results.bins[_mask] ./ _scale)
    _y = vcat(_y, results.jtx[dx, _mask] .* _scale ./ sum(results.jtx[dx, :]))

    # .scaling for diffusive motion
    # _scale = (sqrt(3) * Float64(dx))^(1)

    # _mask = results.jtx[dx, :] .> 0
    # _x = vcat(_x, results.bins[_mask] .^(-3/2))
    # _y = vcat(_y, results.jtx[dx, _mask] .* _scale ./ sum(results.jtx[dx, :]))
  end

  # ?what if i don't bin
  # return _x, _y

  # binEdges = LinRange(extrema(xdata)..., nbins) |> collect
  binEdges = collect(logspace(log.(2, extrema(_x))..., nbins; base=2.0))
  x_binned, y_binned, bin_counts = binning(_x, _y, binEdges)

  _mask = bin_counts .> 0
  return x_binned[_mask], y_binned[_mask]
end
end

if abspath(PROGRAM_FILE) == @__FILE__
  namedtuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

  opts = let
    sts = ArgParseSettings()
    add_arg_table!(sts, "--parameter_space_table", Dict(:arg_type => String, :required => true))
    add_arg_table!(sts, "--nbins", Dict(:arg_type => Int, :default => 4000))
    add_arg_table!(sts, "--cutoff_fraction", Dict(:arg_type => Float64, :default => 0.2))
    add_arg_table!(sts, "--rewrite", Dict(:action => :store_true))
    namedtuple(parse_args(sts))
  end

  df = CSV.read(opts.parameter_space_table, DataFrame)
  @threads for path in df.path
    #. check for existing data, otherwise, if rewrite request, replace data
    if ~ispath(joinpath(path, "processed_ancestry.jld2")) || opts.rewrite
      result::Ancestry.AncestryTools.Results = Ancestry.mainProcess(path, opts.nbins, opts.cutoff_fraction)

      # .save results into database
      jldsave(joinpath(path, "processed_ancestry.jld2"); opts=opts, result=result)
    end
  end
end
