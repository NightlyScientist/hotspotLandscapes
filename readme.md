# Range Expansions across Landscapes with Quenched Noise

This is the public Git repository containing source code used in J. Gonzalez Nuñez, J. Paulose, W. Möbius, D. A. Beller, "Range Expansions across Landscapes with Quenched Noise." arXiv:2310.11563 (2023).

Example usage is shown bellow.

A two-parameter sweep can be performed using the src/processing/generate_parameters_space_search.py file. For example,

> python src/processing/generate_parameter_space_search.py --numberTrials 200 --numberSamples 50 --dims 2000,1100 --data_path workspace/experiments/ --radius 10 --density 0.09 --intensity 8 --ref_line 1000 --detailed_analytics --nEnvs 20 --intervals_1 0.55,0.05,0.8 --intervals_2 0,0.1,1 --rewrite

will generate 20 distinct landscapes and 200 simulations for each landscape, and vary the density (intervals_1) and intensity (intervals_2) for the specified range. A list of options can be found using

> python src/processing/generate_parameter_space_search.py --help

The commands used in the above-referenced article are contained in workarea/simulation_commands.txt. There are two sets of commands, the first corresponds to a parameter scan in the low intensity regime, while the second set of commands correspond to the intermediate to high intensity regime.

Source code used to generate article figures are contained in the figures/ directory. Analysis is performed by providing a list of directories (input variable) containing the simluation data from generate_parameters_space_search.py. If all data saved in workarea/experiments/, the figure scripts will automatically scan and collect all simulation data and proceed with the analysis.

Figures 3A, 7A use data generated from the first set of commands in workarea/simulation_commands.txt, while figures 3B, 7B, 6C, 6D, 9 use both sets of data generation commands. Data for figures 4, 6A, 6B, and 8 are selected from both sets using the parameter values shown in the main text.

Additinally, individual simlalations can be executed using main.jl; a list of command-line options can be displayed using

> julia src/base/main.jl --help
