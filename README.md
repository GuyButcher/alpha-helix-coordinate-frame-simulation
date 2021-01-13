# Readme
### Repository for alpha-helix coordinate frame simulation framework

#### Key Files
* `startup.m` generates any required data `.mat` files and add relevant folders to Matlab's current path.
* `calc_all_bundle_forces.m` runs the full electrostatic interaction simulation and prints the results. This is the exact code used to generate the data used within the manuscript.
* `all_figures.m` When run will generate all figures used within the manuscript along with a supplementary animation. These generated files are save to `methods_paper_figures`.

#### Folders
 * `classes/` holds the class definitions and member functions for the key classes the handle all simulation objects and generate simulation objects from pdb data.
 * `data_generation/` holds the script that generates the object that stores all amino-acid properties used within all parts of the simulation, data generation and figure generatio.
 * `funtions_and_scripts/` holds helper functions.
 * `methods_paper_figures/` holds additional files related to the generation of figures related to the manuscript. It is also the directory that all generated figures will be saved to.
 * `pdb_data/` holds the script and data file that is generated from said script. This script holds to indeces which indicate where each helix starts and ends within their pdb files.
 * `sidechain_dist/` holds the functions that are used to calculate the backbone and sidechain statistical data used to verify our approximations.
 * `talin_pdb_files/` contains the local copy of the pdb files used within the simulation. Each file is originally from SWISS-MODEL and other than being renamed for ease of use, are untouched.
