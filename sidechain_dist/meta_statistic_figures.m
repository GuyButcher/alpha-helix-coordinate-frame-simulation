%% Generate Data

silent_startup;

helix_counts = zeros(13,1);
for i = 1:13
    temp = helix_indeces(:,:,1);
    if (temp(i,5) == 0)
        helix_counts(i) = 4;
    else
        helix_counts(i) = 5;
    end
end
helix_counts_sum = sum(helix_counts);

% Determine number of helices for array of structures pre-alloc

helix_count = 0;
for i = 1:13
    for j = 1:helix_counts(i)
        helix_count = helix_count + 1;
    end
end
% disp(helix_count);

helix_statistics.domain_index = 0;
helix_statistics.helix_index = 0;
helix_statistics.helix_length = 0;
helix_statistics.distances = [];
helix_statistics.fit_line = [];
helix_statistics.distance_mean = 0;
helix_statistics.normalised_distances = [];
helix_statistics.distance_sd = 0;
helix_statistics.distance_variance = 0;
helix_statistics.distance_min = 0;
helix_statistics.distance_max = 0;
helix_statistics.distance_diff = 0;
helix_statistics.fit_line_min = 0;
helix_statistics.fit_line_max = 0;
helix_statistics.fit_line_first = 0;
helix_statistics.fit_line_last = 0;
helix_statistics.fit_line_grad = 0;
helix_statistics.fit_line_angle_d = 0;
helix_statistics_data(helix_count) = helix_statistics;

count = 1;
atom = "CA";
nbins = 10;
progress_bar = waitbar(count / helix_counts_sum, "Calculating Helix Alignment Data");
for i = 1:length(helix_counts)
    domain_index = i;
    for j = 1:helix_counts(i)
        helix_index = j;
        
        [x, distances, fit_line, y_mean_line] = backbone_alignment_analysis(pdb_filepaths,domain_index,helix_indeces,helix_index, atom);
        helix_len = helix_length(pdb_filepaths,domain_index,helix_indeces,helix_index);

        helix_statistics.domain_index = domain_index;
        helix_statistics.helix_index = helix_index;
        helix_statistics.helix_length = helix_len;
        helix_statistics.distances = distances;
        helix_statistics.fit_line = fit_line;
        helix_statistics.distance_mean = y_mean_line(1);
        helix_statistics.normalised_distances = helix_statistics.distances - helix_statistics.distance_mean;
        helix_statistics.distance_sd = std(helix_statistics.normalised_distances);
        helix_statistics.distance_variance = var(helix_statistics.normalised_distances);
        helix_statistics.distance_min = min(helix_statistics.distances);
        helix_statistics.distance_max = max(helix_statistics.distances);
        helix_statistics.distance_diff = helix_statistics.distance_max - helix_statistics.distance_min;
        helix_statistics.fit_line_min = min(helix_statistics.fit_line);
        helix_statistics.fit_line_max = max(helix_statistics.fit_line);
        helix_statistics.fit_line_first = helix_statistics.fit_line(1);
        helix_statistics.fit_line_last = helix_statistics.fit_line(end);
        helix_statistics.fit_line_grad = (helix_statistics.fit_line_last - helix_statistics.fit_line_first) / helix_statistics.helix_length;
        helix_statistics.fit_line_angle_d = atan2d(helix_statistics.fit_line_last - helix_statistics.fit_line_first, helix_statistics.helix_length);

        helix_statistics_data(count) = helix_statistics;
        waitbar(count / helix_counts_sum,progress_bar, sprintf('Calculating Helix Alignment Data: %.0f%%', (count / helix_counts_sum) * 100));
%         fprintf("Count: %d\n", count);
%         fprintf("Index: %d, %d\n", domain_index, helix_index);
        count = count + 1;
    end
end
delete(progress_bar);