function helix_length = helix_length(pdb_filepaths,domain_index,helix_indeces,helix_index)
    chain = Amino_Acid_Chain(pdb_filepaths(domain_index));
    indeces = [helix_indeces(domain_index,helix_index,1);helix_indeces(domain_index,helix_index,2)];
    chain.Select_Data(indeces);
    bb_atoms = Amino_Acid_Chain.Get_Backbone_Atom_Positions(pdb_filepaths(domain_index),indeces(1),indeces(2));
    backbone_positions = chain.backbone_positions(indeces(1):indeces(2),:);
    helix_length = norm(backbone_positions(end,:) - backbone_positions(1,:));
end

