# Range Expansions across Landscapes with Quenched Noise

This is the public Git repository containing source code used in J. Gonzalez Nuñez, J. Paulose, W. Möbius, D. A. Beller, "Range Expansions across Landscapes with Quenched Noise." [arXiv:2310.11563](https://arxiv.org/abs/2310.11563) (2023).

Example usage is shown below.

A two-parameter sweep can be performed using the [src/processing/generate_parameters_space_search.py](./src/processing/generate_parameter_space_search.py) file. For example,

> python src/processing/generate_parameter_space_search.py --numberTrials number_of_independent_runs --numberSamples number_of_sampling_points_for_speed_measure --dims width,height --data_path path_to_save_location --radius radius --density density --intensity intensity --ref_line cutoff_height --detailed_analytics --nEnvs number_of_landscapes --parameters parameter_1,parameter_2 --intervals_1 start,step,stop --intervals_2 start,step,stop --rewrite

will generate `number_of_landscapes` distinct landscapes and `number_of_trials` simulations for each landscape, and vary the `parameter_1` and `parameter_2` parameters for the specified range (`interval_1` and `interval_2`); data will be saved to `path_to_save_loation`; hotspot parameters are set to `density`, `intensity`, and `radius` (see below); system size is set to `width` x `height`; the simulation is terminated at `cutoff_height`; measurements of front propgation speed require setting the number of sampling points during a simulation to `number_of_sampling_points_for_speed_measure`. 

Note that `parameter_1` and `parameter_2` must be one of the following strings: ["density", "intensity", "radius"], for example, `--parameters "density","intensity"`. While values for all three `--density`, `--intensity`, `--radius` can be provided, only the parameter that will not be varied is required, for example "radius"; the other two values, if provided, will be overwritten by values set by `interval_1` and `interval_2`.

A list of simulation options and their descriptions can be found using

> python src/processing/generate_parameter_space_search.py --help

The commands used in the above-referenced article are contained in [workarea/simulation_commands.txt](./workarea/simulation_commands.txt). There are two sets of commands; the first corresponds to a parameter scan in the low intensity regime, while the second set of commands corresponds to the intermediate to high intensity regime.

Source code used to generate article figures are contained in the [figures/](./figures/) directory. Analysis is performed by providing a list of directories (input variable) containing the simluation data from [generate_parameter_space_search.py](./src/processing/generate_parameter_space_search.py). If all data are saved in, for example, workarea/experiments/, the figure scripts will automatically scan and collect all simulation data and proceed with the analysis.

There are three pre-processing scripts that need to be executed before running figure scripts: [figures/routines/ancestry.jl](figures/routines/ancestry.jl) which builds a n-ary tree using lineage branch points and processes lineage coalescences; [src/processing/process_fastest_paths.jl](src/processing/process_fastPaths.jl) which constructs a hotspot graph and uses a Floyd-Warshall (Dijkstra) path-finding algorithm to find optimal paths; [src/processing/process_lineages.jl](src/processing/process_lineages.jl) which measures the number of surviving ancestors, lineage Mean-Square-Displacement, and lineage tortuosity. The command-line options for these scripts can be listed using the `--help` option. 

Figures 3A and 7A in the article use data generated from the first set of commands in [workarea/simulation_commands.txt](./workarea/simulation_commands.txt), while figures 3B, 7B, 6C, 6D, and 9 use both sets of data generation commands. Data for figures 4, 6A, 6B, and 8 are selected from both sets using the parameter values shown in the main text.

Additionally, individual simulations can be executed using main.jl; a list of command-line options can be displayed using

> julia src/base/main.jl --help

SI figures of optimal path calculations can be generated in the jupyter notebook [src/calculations/hotspot-graph.ipynb](./src/calculations/hotspot-graph.ipynb). 

This code requires both Julia (v1.10.4) and python (v3.11.5) to be installed. Our Julia environment is contained in the [Manifest.toml](Manifest.toml) and [Project.toml](Project.toml) files. Our Python environment is provided in [environment.yml](environment.yml).

While the script [generate_parameter_space_search.py](./src/processing/generate_parameter_space_search.py) to generate simulations works best with the [slurm workload manager](https://slurm.schedmd.com/overview.html) installed, the script will check for an existing slurm installion and will fallback to executing sequentially via bash if no slurm installation is found.
