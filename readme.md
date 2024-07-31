# Range Expansions across Landscapes with Quenched Noise

This is the public Git repository containing source code used in J. Gonzalez Nuñez, J. Paulose, W. Möbius, D. A. Beller, "Range Expansions across Landscapes with Quenched Noise." arXiv:2310.11563 (2023).

Example usage is shown bellow.

A two-parameter sweep can be performed using the src/processing/generate_parameters_space_search.py file. For example,

> python src/processing/generate_parameter_space_search.py --numberTrials number_of_independent_runs --numberSamples number_of_sampling_points_for_speed_measure --dims width,height --data_path path_to_save_location --radius radius --density density --intensity intensity --ref_line cutoff_height --detailed_analytics --nEnvs number_of_landscapes --intervals_1 start,step,stop --intervals_2 start,step,stop --rewrite

will generate number_of_landscapes distinct landscapes and number_of_trials simulations for each landscape, and vary the density (intervals_1) and intensity (intervals_2) for the specified range; data will be saved to path_to_save_loation; hotspot parameters are set to density, intensity radius; system size is set to width x height; simulation data is terminated at cutoff_height; measurements of front propgation speed require setting a the number of sampling points during a simulation number_of_sampling_points_for_speed_measure. 

A list of simulation options can be found using

> python src/processing/generate_parameter_space_search.py --help

The commands used in the above-referenced article are contained in workarea/simulation_commands.txt. There are two sets of commands, the first corresponds to a parameter scan in the low intensity regime, while the second set of commands correspond to the intermediate to high intensity regime.

Source code used to generate article figures are contained in the figures/ directory. Analysis is performed by providing a list of directories (input variable) containing the simluation data from generate_parameters_space_search.py. If all data saved in, for example, workarea/experiments/, the figure scripts will automatically scan and collect all simulation data and proceed with the analysis.

Figures 3A, 7A use data generated from the first set of commands in workarea/simulation_commands.txt, while figures 3B, 7B, 6C, 6D, 9 use both sets of data generation commands. Data for figures 4, 6A, 6B, and 8 are selected from both sets using the parameter values shown in the main text.

Additinally, individual simlalations can be executed using main.jl; a list of command-line options can be displayed using

> julia src/base/main.jl --help
