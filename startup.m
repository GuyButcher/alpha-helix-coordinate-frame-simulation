clc
fprintf("---: PDB and Talin Sim :---\n")
fprintf("Executing Startup:...\n");

fprintf("Adding 'classes' path.\n");
addpath("classes/");

fprintf("Adding 'data_generation' path.\n");
addpath("data_generation/");

fprintf("Running 'gen_residue_properties'.\n");
gen_residue_properties;

fprintf("Adding 'functions_and_scripts' path.\n");
addpath("functions_and_scripts/");

fprintf("Adding 'pdb_data' path.\n");
addpath("pdb_data/");

fprintf("Adding 'methods_paper_figures' path.\n");
addpath("methods_paper_figures//");

fprintf("Adding 'sidechain_dist' path.\n");
addpath("sidechain_dist//");

fprintf("Loading variable 'helix_indeces' to Workspace.\n");
if(~isfile('pdb_data\helix_indeces.mat'))
    generate_data;
end
load("helix_indeces.mat");

fprintf("Loading variable 'pdb_filepaths' to Workspace.\n")
pdb_filepaths = [ ...
    "talin_pdb_files/R1.pdb",...
    "talin_pdb_files/R2.pdb",...
    "talin_pdb_files/R3.pdb",...
    "talin_pdb_files/R4.pdb",...
    "talin_pdb_files/R5.pdb",...
    "talin_pdb_files/R6.pdb",...
    "talin_pdb_files/R7R8.pdb",...
    "talin_pdb_files/R7R8.pdb",...
    "talin_pdb_files/R9.pdb",...
    "talin_pdb_files/R10.pdb",...
    "talin_pdb_files/R11.pdb",...
    "talin_pdb_files/R12.pdb",...
    "talin_pdb_files/R13.pdb",...
    ];


fprintf("Startup Complete.\n\n\n");