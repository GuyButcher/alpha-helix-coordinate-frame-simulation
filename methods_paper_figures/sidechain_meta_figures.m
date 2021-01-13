%% Scatter Graph of each type of amino acid with ca to sidechain centre

data_structs = sidechain_atom_extraction();

amino_acid_index = 2;

num_domains = length(data_structs);
num_amino_acids = length(data_structs{1});

transformed_sidechain_atom_positions = cell(num_amino_acids,num_domains);
transformed_sidechain_atom_distances = cell(num_amino_acids,num_domains);

for domain_index = 1:num_domains

    current_residue = data_structs{domain_index}{amino_acid_index};

    num_residues = length(current_residue);
    num_atoms = length(current_residue(1).raw_sidechain_positions(:,1));

    transformed_atom_positions = zeros(num_atoms,3,num_residues);
    transformed_atom_distances = zeros(num_atoms,num_residues);

    for residue_index = 1:num_residues
        ca_position = current_residue(residue_index).ca_position;
        raw_position = current_residue(residue_index).raw_sidechain_positions;
        mean_raw_position = mean(raw_position);
        normalised_position = raw_position - ca_position;
        mean_normalised_position = mean(normalised_position);
        vector_ca_to_mean_position = (mean_raw_position - ca_position) / vecnorm(mean_raw_position - ca_position);
        vector_ca_to_mean_position = (mean_normalised_position) / vecnorm(mean_normalised_position);
        vector_x_axis = [1 0 0];
        rot = VecToVecRotation(vector_ca_to_mean_position, vector_x_axis);
        transformed_atom_position = (rot * normalised_position')';
%         transformed_atom_position = transformed_atom_position - ca_position;
%         transformed_atom_position = transformed_atom_position - mean(transformed_atom_position);
        
        transformed_atom_positions(:,:,residue_index) = transformed_atom_position;
        transformed_atom_distances(:,residue_index) = vecnorm(transformed_atom_position')';

    end
    
    transformed_sidechain_atom_positions{amino_acid_index,domain_index} = transformed_atom_positions;
    transformed_sidechain_atom_distances{amino_acid_index,domain_index} = transformed_atom_distances;

end

%%
figure();
hold on
for domain_index = 1:num_domains
    num_residues = length(transformed_sidechain_atom_positions{amino_acid_index,domain_index}(1,1,:));
    for residue_index = 1:num_residues
        scatter3( ...
            transformed_sidechain_atom_positions{amino_acid_index,domain_index}(:,1,residue_index), ...
            transformed_sidechain_atom_positions{amino_acid_index,domain_index}(:,2,residue_index), ...
            transformed_sidechain_atom_positions{amino_acid_index,domain_index}(:,3,residue_index),'.');
    end
end
xlabel("X");
ylabel("Y");
zlabel("Z");
hold off
%%
figure();
atom_distances = [];
for domain_index = 1:num_domains
    current_distance = transformed_sidechain_atom_distances{amino_acid_index,domain_index}(:);
    atom_distances = [atom_distances;current_distance];
end
histfit(atom_distances);

%% Data Generation
data_structs = sidechain_atom_extraction();

%%
num_domains = length(data_structs);
num_amino_acids = length(data_structs{1});

sidechain_atom_positions = cell(num_amino_acids,num_domains);
sidechain_atom_distances = cell(num_amino_acids,num_domains);

for amino_index = 1:num_amino_acids

    for domain_index = 1:num_domains

        current_residue = data_structs{domain_index}{amino_index};
        
        if (isempty(current_residue)); continue; end

        num_residues = length(current_residue);
        num_atoms = length(current_residue(1).raw_sidechain_positions(:,1));

        atom_positions = zeros(num_atoms,3,num_residues);
        atom_distances = zeros(num_atoms,num_residues);

        for residue_index = 1:num_residues
            ca_position = current_residue(residue_index).ca_position;
            raw_position = current_residue(residue_index).raw_sidechain_positions;
            ca_corrected_position = raw_position - ca_position;
            zero_meaned_position = ca_corrected_position - mean(ca_corrected_position);
            distances = vecnorm(zero_meaned_position')';


            atom_positions(:,:,residue_index) = zero_meaned_position;
            atom_distances(:,residue_index) = distances;
        end

        sidechain_atom_positions{amino_index,domain_index} = atom_positions;
        sidechain_atom_distances{amino_index,domain_index} = atom_distances;


    end
end
%% Figure Generation

% Figure For Arg

% Final Data Needed:
% - Scatter Graph of each type of amino acid without regard for direction
% - Scatter Graph of each type of amino acid with ca to sidechain centre
% aligned
% - Distribution plot of distances from sidechain atom pos to sidechain
% centre pos

% Arg test plot

data_structs = sidechain_atom_extraction();

amino_acid_index = 1;

num_domains = length(data_structs);
num_amino_acids = length(data_structs{1});

sidechain_atom_positions = cell(num_amino_acids,num_domains);
sidechain_atom_distances = cell(num_amino_acids,num_domains);

for domain_index = 1:num_domains

    current_residue = data_structs{domain_index}{amino_acid_index};

    num_residues = length(current_residue);
    num_atoms = length(current_residue(1).raw_sidechain_positions(:,1));

    atom_positions = zeros(num_atoms,3,num_residues);
    atom_distances = zeros(num_atoms,num_residues);

    for residue_index = 1:num_residues
        ca_position = current_residue(residue_index).ca_position;
        raw_position = current_residue(residue_index).raw_sidechain_positions;
        ca_corrected_position = raw_position - ca_position;
        zero_meaned_position = ca_corrected_position - mean(ca_corrected_position);
        distances = vecnorm(zero_meaned_position')';


        atom_positions(:,:,residue_index) = zero_meaned_position;
        atom_positions(:,:,residue_index) = ca_corrected_position;
        atom_distances(:,residue_index) = distances;
    end
    
    sidechain_atom_positions{amino_acid_index,domain_index} = atom_positions;
    sidechain_atom_distances{amino_acid_index,domain_index} = atom_distances;
    

end

% Plot Scatter the ARG sidechain atom positions relative to sidechain
% mid-position.

figure();
hold on
instance_counter = zeros(num_domains,1);
for domain_index = 1:num_domains
    num_residues = length(sidechain_atom_positions{amino_acid_index,domain_index}(1,1,:));
    for residue_index = 1:num_residues
        scatter3( ...
            sidechain_atom_positions{amino_acid_index,domain_index}(:,1,residue_index), ...
            sidechain_atom_positions{amino_acid_index,domain_index}(:,2,residue_index), ...
            sidechain_atom_positions{amino_acid_index,domain_index}(:,3,residue_index),'.');
        instance_counter(domain_index) = instance_counter(domain_index) + 1;
    end
end
hold off


% Plot Scatter the Arg sidechain mid-position relative to CA positions.

figure();
hold on
instance_counter = zeros(num_domains,1);
for domain_index = 1:num_domains
    num_residues = length(sidechain_atom_positions{amino_acid_index,domain_index}(1,1,:));
    for residue_index = 1:num_residues
        scatter3( ...
            sidechain_atom_positions{amino_acid_index,domain_index}(:,1,residue_index), ...
            sidechain_atom_positions{amino_acid_index,domain_index}(:,2,residue_index), ...
            sidechain_atom_positions{amino_acid_index,domain_index}(:,3,residue_index),'.');
        instance_counter(domain_index) = instance_counter(domain_index) + 1;
    end
end
hold off


%%

% Extract distances into 2D array from all domains

amino_acid_distances = [];
for domain_index = 1:num_domains
    current_distances = sidechain_atom_distances{amino_acid_index,domain_index}(:);
    amino_acid_distances = [amino_acid_distances;current_distances];
end
histfit(amino_acid_distances);
boxplot(amino_acid_distances);
mean(amino_acid_distances);

%%

