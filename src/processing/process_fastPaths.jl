using FileIO, JLD2, StatsBase, ArgParse, Base.Threads, DataFrames, CSV
include("../calculations/lineages.jl")
include("../calculations/OptimalPaths.jl")
include("../../figures/common/theme.jl")
include("../../figures/common/modeling.jl")
using .ContinuousOptPaths, DelaunayTriangulation

# namedtuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

# >visualization of voronoi tessellation with nearest neighbors
function voronoiGraphImage(voronoiGraph, opts, save_name)
  fig = Figure(; size=(800, 800), figure_padding=30)
  ax = Axis(fig[1, 1]; title="Voronoi tessellation with graph", titlealign=:left)

  areas = get_area.(Ref(voronoiGraph.voronoi), 1:num_polygons(voronoiGraph.voronoi))
  colors = get(cgrad(ColorSchemes.grayC), areas, :extrema)

  voronoiplot!(ax, voronoiGraph.voronoi; show_generators=true, markercolor=:green, color=colors)

  for edge_candidate in voronoiGraph.edgeGenerators
    scatter!(ax, get_generator(voronoiGraph.voronoi, edge_candidate)...; color=:lightgreen, markersize=15)
  end

  for (gen_index, nbors) in voronoiGraph.connections
    for nbor in voronoiGraph.connections[gen_index]
      x1, y1 = get_generator(voronoiGraph.voronoi, gen_index)
      x2, y2 = get_generator(voronoiGraph.voronoi, nbor)
      abs(x2 - x1) > fld(opts.width, 2) && continue
      lines!(ax, [x1, x2], [y1, y2]; color=:dodgerblue)
    end
  end

  fig_path = save_name * "/figure_voronoiGraph.png"
  limits!(ax, 0, opts.width, 0, opts.ref_line)
  # @info fig_path
  save(fig_path, fig)
end

# > snapshot of optimal paths and lineage traces
function optimalPathsImage(
  sourceIndices, optimalPathSets, voronoiGraph, opts, save_name; vplot=false, lplot=false, pplot=false
)
  input = save_name
  fig = Figure(; size=(800, 600), figure_padding=20)
  ax = Axis(fig[1, 1]; backgroundcolor=:white)

  # .draw voronoi tessellation
  areas = get_area.(Ref(voronoiGraph.voronoi), 1:num_polygons(voronoiGraph.voronoi))
  colors = get(cgrad(ColorSchemes.grayC), areas, :extrema)

  vplot && voronoiplot!(ax, voronoiGraph.voronoi; show_generators=true, markercolor=:purple, color=colors)

  if lplot
    # .draw lineageMap by smoothing lineages with a 7x7 averaging kernel
    results::LineageTracing.Results = LineageTracing.lineageTraces(input, opts; full=false)
    smoothed_lineageMap = copy(results.lineageMap)
    for x in axes(smoothed_lineageMap, 1), y in axes(smoothed_lineageMap, 2)
      _xrange = clamp.((x - 3):(x + 3), 1, opts.width)
      _yrange = clamp.((y - 3):(y + 3), 1, opts.ref_line)
      smoothed_lineageMap[x, y] = sum(results.lineageMap[_xrange, _yrange]) / 49
    end

    smoothed_lineageMap ./= mean(smoothed_lineageMap[smoothed_lineageMap .> 0])
    heatmap!(ax, smoothed_lineageMap; lowclip=:transparent, colorrange=(0.01, 0.5), colormap=ColorSchemes.algae)
  end

  # .draw line segments between hotspots
  distances = Float64[]
  for source_index in sourceIndices
    push!(distances, getfield.(optimalPathSets[source_index], :pathDistance)...)
  end

  cmap = alphaColor(ColorSchemes.berlin, 0.8)
  colors = get(cgrad(cmap), distances, :extrema)

  if pplot
    for source_index in sourceIndices
      for optimalPath in optimalPathSets[source_index]
        color = popfirst!(colors)

        for j in 2:2:length(optimalPath.pathPositions)
          j == firstindex(optimalPath.pathPositions) && continue
          (x1, y1) = optimalPath.pathPositions[j - 1]
          (x2, y2) = optimalPath.pathPositions[j]

          lines!(ax, [x1, x2], [y1, y2]; color=color, linewidth=5)
          # lines!(ax, [x1, x2], [y1, y2]; color=:blue, linewidth=5)
        end
      end
    end
  end

  # for _visited_generatorIndex in _mask
  #   scatter!(ax, generators[_visited_generatorIndex]...; markersize=30, color=:pink)
  # end

  # .draw segments from top edges to hotspots
  # for i in axes(sources_to_htspts, 1)
  #   xposition, yposition = (Int64(i), opts.ref_line)
  #   for j in axes(sources_to_htspts, 2)
  #     centerX, centerY = edgeGencenters[sources_to_htspts[i, j]]
  #     lines!(ax, Int[xposition, centerX], Int[yposition, centerY]; color=(:blue, 0.2))
  #   end
  # end

  # # .draw segments from bottom edges to hotspots
  # for i in axes(sources_to_htspts, 1)
  #   xposition, yposition = (Int64(i), opts.ref_line)
  #   for j in axes(sources_to_htspts, 2)
  #     centerX, centerY = edgeGencenters[sources_to_htspts[i, j]]
  #     lines!(ax, Int[xposition, centerX], Int[yposition, centerY]; color=(:red, 0.2))
  #   end
  # end

  if pplot
    kwargs = (ticksize=15, tickalign=1, ticklabelsize=25, colormap=ColorSchemes.berlin, labelsize=25)
    # ticks = ([0, 0.5, 1], ["0", "0.5", "1"])
    colorrange = extrema(distances)
    Colorbar(fig[1, 2]; label="Path Travel Time", kwargs..., colorrange=colorrange)
  end

  if lplot
    kwargs = (
      ticksize=15,
      tickalign=1,
      ticklabelsize=25,
      colormap=ColorSchemes.algae,
      vertical=false,
      flipaxis=false,
      labelsize=25
    )
    # ticks = ([0, 0.5, 1], ["0", "0.5", "1"])
    colorrange = (0.01, 0.5)
    Colorbar(fig[2, 1]; label="Lineage Visitation Frequency", kwargs..., colorrange=colorrange)
  end

  scalebar!(ax, (110, 100); width=200, height=20, color=:red, fontsize=35)
  limits!(ax, 1, opts.width, 1, opts.ref_line)
  colgap!(fig.layout, 1)
  rowgap!(fig.layout, 1)
  hidedecorations!(ax)
  fig_path = save_name * "/figure_optimalPaths"
  vplot && (fig_path *= "_voronoi")
  lplot && (fig_path *= "_lineages")
  pplot && (fig_path *= "_optimals")
  save(fig_path * ".png", fig)
  # @info fig_path
end

function execute_main(input, n_fastest, cutoff, interval; NN = false, figures=false)
  save_name = input
  opts = namedtuple(load(joinpath(input, "Opts.jld2")))
  htspts = load(joinpath(input, "objects.jld2"), "objs")

  # .voronoi triangulation
  htspts = filter(t -> last(t) <= opts.ref_line, unique!(htspts))
  points = transpose(hcat(Float64.(first.(htspts)), Float64.(last.(htspts))))

  voronoiGraph = ContinuousOptPaths.voronoiTriangulation(points, opts)
  figures && voronoiGraphImage(voronoiGraph, opts, save_name)

  # .nearest hotspot to each source (sink) node
  sinks_to_htpts, sources_to_htspts = ContinuousOptPaths.sourceSinkConnections(
    voronoiGraph, opts; ref_line=opts.ref_line, n_near=30, interval=interval
  )

  # .combine sources, sinks, and hotspots as nodes in the graph
  g = weightedGraph(voronoiGraph, sources_to_htspts, sinks_to_htpts, opts; NN=NN)

  sourceIndices, sinkIndices = ContinuousOptPaths.edgeNodeIndices(
    sources_to_htspts, sinks_to_htpts, voronoiGraph.generators
  )
  gPositions = ContinuousOptPaths.graphPosisions(sources_to_htspts, sinks_to_htpts, voronoiGraph.generators, opts)

  floyd = NN ? false : true
  optimalPathSets = allOptimalPaths(
    sourceIndices, sinkIndices, opts, gPositions, g; n=n_fastest, cutoff=cutoff, floyd=floyd
  )

  if figures
    optimalPathsImage(sourceIndices, optimalPathSets, voronoiGraph, opts, save_name; pplot=true, vplot=true)
    optimalPathsImage(sourceIndices, optimalPathSets, voronoiGraph, opts, save_name; lplot=true, pplot=true)
    optimalPathsImage(sourceIndices, optimalPathSets, voronoiGraph, opts, save_name; lplot=true, vplot=true)
  end

  #.visited hotspots, MSD, and suriving ancestors
  visitedObjectsIndices = ContinuousOptPaths.objectVisits(voronoiGraph, sourceIndices, optimalPathSets)
  opt_results = ContinuousOptPaths.optimalPathMSD(optimalPathSets, opts.ref_line, opts.width)
  survivingAncestors = ContinuousOptPaths.pathRoots(optimalPathSets)

  # .save data to csv: distance | mean 
  CSV.write(
    joinpath(input, "processed_optimalPathMSD.csv"),
    DataFrame(; distance=collect(1:(opts.ref_line)), msd=opt_results.MSD)
  )
  jldsave(joinpath(input, "processed_visitedObjects_optimalPaths.jld2"); visitedObjectsIndices=visitedObjectsIndices)
  jldsave(joinpath(input, "processed_surivingAncestors_optimalPaths.jld2"); survivingAncestors=survivingAncestors)
end

opts = begin
  sts = ArgParseSettings()
  add_arg_table!(sts, "--parameter_space_table", Dict(:arg_type => String, :required => true))
  add_arg_table!(sts, "--n_fastest", Dict(:arg_type => Int, :required => true))
  add_arg_table!(sts, "--interval", Dict(:arg_type => Int, :required => true))
  add_arg_table!(sts, "--cutoff", Dict(:arg_type => Float64, :required => true))
  add_arg_table!(sts, "--figures", Dict(:action => :store_true))
  add_arg_table!(sts, "--NN", Dict(:action => :store_true))
  parse_args(sts)
end

df = CSV.read(opts["parameter_space_table"], DataFrame)

# task: option to process only a small slice through the phase space

@threads for single_trial in eachrow(df)
  data_path = single_trial.path
  execute_main(data_path, opts["n_fastest"], opts["cutoff"], opts["interval"]; figures=opts["figures"], NN=opts["NN"])
end
