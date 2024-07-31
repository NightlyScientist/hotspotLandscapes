using DataFrames, CSV, ArgParse, JLD2, FileIO, Base.Threads, StatsBase
include("../../src/calculations/lineages.jl")

function execute_main(window_method, data_path, height, width, ref_line, num_trials; interval=1)
  method = window_method ? LineageTracing.lineageMSD_window : LineageTracing.lineageMSD
  y_max = window_method ? 750 : height

  ancestorCounts = zeros(Int64, width, num_trials)
  msd = zeros(Float64, (y_max, num_trials))
  tortuosity = zeros(Float64, num_trials)

  _processed = zeros(Bool, num_trials)
  sources = LineageTracing.graphSources(ref_line, width, height)

  dataFile = joinpath(data_path, "data_phylo.jld2")
  isfile(dataFile) || return nothing

  jldopen(dataFile, "r") do file
    trials = LineageTracing.mask(keys(file), "trial")
    for j in eachindex(trials)
      phylo = file[trials[j]]["phylogeny"]
      results = method(height, width, ref_line, sources, phylo; interval=interval)
      msd[:, j] .= results.MSD
      tortuosity[j] = first(results.tortuosity)

      _uniqueLabels = LineageTracing.phylogenyTreeRoots(phylo, width)
      ancestorCounts[_uniqueLabels, j] .+= 1
      _processed[j] = true
    end
  end

  index_mask = findall(_processed)
  isempty(index_mask) && return nothing

  # .save tortuosity
  write(joinpath(data_path, "processed_lineage_tortuosity.txt"), join(mean_and_std(tortuosity[index_mask]), ","))

  # .save data to csv: matrix -> dataframe: cols = trials, rows = x position
  CSV.write(
    joinpath(data_path, "processed_lineage_ancestorCounts.csv"), DataFrame(ancestorCounts[:, index_mask], :auto)
  )

  # .save data to csv: mean | std | stderr
  CSV.write(
    joinpath(data_path, "processed_lineageMSD.csv"),
    DataFrame(;
      distance=collect(1:y_max),
      msd=mean(msd[:, index_mask]; dims=2)[1:y_max],
      std_error=std(msd[:, index_mask]; dims=2)[1:y_max] ./ sqrt(num_trials)
    )
  )
end

opts = let
  sts = ArgParseSettings()
  add_arg_table!(sts, "--parameter_space_table", Dict(:arg_type => String, :required => true))
  add_arg_table!(sts, "--window_method", Dict(:action => :store_true))
  add_arg_table!(sts, "--interval", Dict(:arg_type => Int, :default=>1))
  parse_args(sts)
end

df = CSV.read(opts["parameter_space_table"], DataFrame)

width = df[1, :].width
height = df[1, :].height
dims = (width, height)

@threads for single_trial in eachrow(df)
  ref_line = single_trial.ref_line
  num_trials = single_trial.numberTrials
  data_path = single_trial.path
  execute_main(opts["window_method"], data_path, height, width, ref_line, num_trials; interval=opts["interval"])
end
