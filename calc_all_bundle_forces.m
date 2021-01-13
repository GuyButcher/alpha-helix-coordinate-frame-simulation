clc
clear
startup

fprintf('Running Sim:\n');
n = 1;
for n = 1:13
    pdb_path = pdb_filepaths(n);
    indeces = helix_indeces(n,:,:);
    tic
    sim = Simulation();
    [sum, forces, combi] = sim.Run_Bundle_Static_Force(pdb_path,indeces);
    toc
    output(n).sum = sum;
    output(n).forces = forces;
    output(n).combination = combi;
end
fprintf('Done!\n');
[output.sum]'