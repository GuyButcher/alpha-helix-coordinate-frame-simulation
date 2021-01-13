function sc_dist(pdb_file_index, amino_acid_index, nbins)
    silent_startup;
    names = load('sidechain_dist/amino_acid_names.mat');
    names = names.names;
    
    include_mean_dist = true;
    
    colourful = false;
    
    if (nargin < 3)
        nbins = 10;
    end
    
    
    num_aa_t = length(names);
    aa_positions = cell(num_aa_t,1);
    bb_positions = cell(num_aa_t,1);

    filename = convertStringsToChars(pdb_filepaths(pdb_file_index));

    pdb_file = pdbread(filename);
    atoms = pdb_file.Model.Atom;
    for i = 1:length(atoms)
        atoms(i).AtomName = convertCharsToStrings(atoms(i).AtomName);
    end


    residue_indeces = unique([atoms.resSeq]);
    residue_count = length(residue_indeces);
    aminoNames = strings(residue_count,1);
    aminos = zeros(residue_count,3);
    sidechains = zeros(residue_count,3);
    sidechains_sizes = zeros(residue_count,1);
    backbone_positions = zeros(residue_count,3);
    backbone_position_distances = zeros(residue_count,1);

    for i = 1:residue_count
        index = [atoms.resSeq] == residue_indeces(i);
        sequence = atoms(index);

        aminoNames(i) = convertCharsToStrings(sequence(2).resName);
        aminos(i,1) = sequence(2).X;
        aminos(i,2) = sequence(2).Y;
        aminos(i,3) = sequence(2).Z;
        backbone_positions(i,:) = aminos(i,:);
        r = find([sequence.AtomName] == "CB");
        if (r ~= 0)

            sequence([sequence.AtomName] == "N") = [];
            sequence([sequence.AtomName] == "CA") = [];
            sequence([sequence.AtomName] == "C") = [];
            sequence([sequence.AtomName] == "O") = [];
            num_atoms = length(sequence);

            if(any([sequence.AtomName] == "OXT"))
                sequence([sequence.AtomName] == "OXT") = [];
                %fprintf("Removing N-Terminus Atom.\n");
            end

            sidechains(i,1) = sum([sequence.X])/num_atoms;
            sidechains(i,2) = sum([sequence.Y])/num_atoms;
            sidechains(i,3) = sum([sequence.Z])/num_atoms;

            aa_index = (names == aminoNames(i));
            position = cat(2,[sequence.X]',[sequence.Y]',[sequence.Z]') - sidechains(i,:);
            aa_positions{aa_index} = cat(3,aa_positions{aa_index},position);
            bb_positions{aa_index} = cat(3,bb_positions{aa_index},aminos(i,:) - sidechains(i,:));

    %         distances = zeros(num_atoms,1);
    % 
    %          for j = 1:num_atoms
    %             distances(j) = sqrt( ... 
    %                 (sidechains(i,1) - sequence(j).X)^2 + ...
    %                 (sidechains(i,2) - sequence(j).Y)^2 + ...
    %                 (sidechains(i,3) - sequence(j).Z)^2);
    %          end      
        end
    end
    
    
    

    aa_count = length(aa_positions{amino_acid_index}(1,1,:));
    aa_atom_count = length(aa_positions{amino_acid_index}(:,1,1));
    aa_x_cat = [];
    aa_y_cat = [];
    aa_z_cat = [];
    
   
    bb_distances = [];
    for i = 1:aa_count
        mean_sc_pos = mean(aa_positions{amino_acid_index}(:,:,i));
        bb_distances(i) = norm(bb_positions{amino_acid_index}(:,:,i) - mean_sc_pos);
    end

    aa_distances = [];
    for i = 1:aa_count
        for j = 1:aa_atom_count
            aa_distances = cat(1,aa_distances,norm(aa_positions{amino_acid_index}(j,:,i)));
        end
    end
    for i = 1:aa_count
        aa_x_cat = cat(1,aa_x_cat,aa_positions{amino_acid_index}(:,1,i));
        aa_y_cat = cat(1,aa_y_cat,aa_positions{amino_acid_index}(:,2,i));
        aa_z_cat = cat(1,aa_z_cat,aa_positions{amino_acid_index}(:,3,i));
    end
    figure();
    hold on
    if ~colourful
    plot3(aa_x_cat,aa_y_cat,aa_z_cat,'rx');
    else
        for i = 1:length(aa_positions{amino_acid_index}(1,1,:))
            plot3(aa_positions{amino_acid_index}(:,1,i),...
                aa_positions{amino_acid_index}(:,2,i),...
                aa_positions{amino_acid_index}(:,3,i),'x','Color',rand(1,3));
        end
    end
    plot3(0,0,0,'og');
    limitsy = ylim;
    if(include_mean_dist == true)
        label = plot3([mean(bb_distances),mean(bb_distances)],[limitsy(1),limitsy(2)],[0,0]);
        legend(label,'Mean distance to Backbone');
    end
    title(strcat("R",num2str(pdb_file_index), ", ", convertCharsToStrings(names(amino_acid_index))));
    axis equal
    hold off
    figure();
    % hist = histogram(aa_distances,25);
    hist = histfit(aa_distances,nbins);
    title(strcat("R",num2str(pdb_file_index), ", ", convertCharsToStrings(names(amino_acid_index))));
end