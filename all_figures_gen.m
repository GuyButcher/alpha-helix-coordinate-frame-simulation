%% 3D Helix visualisations

startup;

sim = Load_Sim_Talin_Subdomain(8);
fig = sim.Make_Figure_Dots();
ax = gca();
view(ax,-12,30);
lgt = camlight('Headlight');
drawnow();
print('methods_paper_figures/figure_a_subdomain_bundle_view','-dpng','-r300');

%%

startup;

sim = Load_Sim_Talin_Subdomain(8);
fig = sim.Make_Figure();
ax = gca();
view(ax,-12,30);
lgt = camlight('Headlight');
drawnow();
print('methods_paper_figures/figure_b_subdomain_bundle_view','-dpng','-r300');

%% Backbone approximation scatter and histfit plots

close all

bp_atoms_dist(7,2,8);

fig1 = figure(1);
ax1 = fig1.CurrentAxes;
title(ax1,'');


fig2 = figure(2);
ax2 = fig2.CurrentAxes;
title(ax2,'');

drawnow

print(fig1, 'methods_paper_figures\figure_bb_scatter_1','-dpng','-r300');
print(fig2, 'methods_paper_figures\figure_bb_hist_1','-dpng','-r300');

%% Backbone approximation scatter and histfit plots
close all

bp_atoms_dist(4,3,8);

fig1 = figure(1);
ax1 = fig1.CurrentAxes;
title(ax1,'');


fig2 = figure(2);
ax2 = fig2.CurrentAxes;
title(ax2,'');

drawnow

print(fig1, 'methods_paper_figures\figure_bb_scatter_2','-dpng','-r300');
print(fig2, 'methods_paper_figures\figure_bb_hist_2','-dpng','-r300');

%% Backbone meta statistics plots

close all

meta_statistic_figures;

fitted_angles = [helix_statistics_data.fit_line_angle_d]';
% disp(horzcat((1:61)',fitted_angles));
fitted_angles_mean = mean(fitted_angles);
fitted_angles_sd = std(fitted_angles - fitted_angles_mean);
fitted_angles_min = min(fitted_angles);
fitted_angles_max = max(fitted_angles);
fprintf("Normalised Angles of Fitted line: ");
fprintf("Min: %f\nMax: %f\nMean: %f\nStandard Deviation: %f\n", fitted_angles_min, fitted_angles_max, fitted_angles_mean, fitted_angles_sd);
boxplot(fitted_angles, 'Labels', "All");
ylabel("Angle of fitted line, Degrees");
title("Boxplot - Angles of line fitted to backbone position distribution");

drawnow

fig1 = figure(1);
print(fig1, 'methods_paper_figures\backbone-skew-boxplot', '-dpng', '-r150');

figure();
histfit(fitted_angles, 15);
title('Distribution of fitted line angles');
ylabel('Binned Count');
xlabel('Angle of Regression line, Degrees');

drawnow

fig2 = figure(2);
print(fig2, 'methods_paper_figures\backbone-skew-histogram', '-dpng', '-r150');

%% Sidechain approximaation scatter and histfit a_i 10

close all

amino_index = 10;

data_structs = sidechain_atom_extraction();
num_domains = length(data_structs);
input_data = process_sidechain_data(data_structs);

az = 0;
el = 0;

output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);
fig1_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax1 = fig1_handle.CurrentAxes;
set(ax1, 'DataAspectRatio', [1 1 1]);
view(ax1,az,el);
title(ax1,'0 Degrees Rotation');

az = 30;
el = 0;

fig2_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax2 = fig2_handle.CurrentAxes;
set(ax2, 'DataAspectRatio', [1 1 1]);
view(ax2,az,el);
title(ax2,'30 Degrees Rotation');


az = 60;
el = 0;

fig3_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax3 = fig3_handle.CurrentAxes;
set(ax3, 'DataAspectRatio', [1 1 1]);
view(ax3,az,el);
title(ax3,'60 Degrees Rotation');

az = 90;
el = 0;

fig4_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax4 = fig4_handle.CurrentAxes;
set(ax4, 'DataAspectRatio', [1 1 1]);
view(ax4,az,el);
title(ax4,'90 Degrees Rotation');

drawnow
% 
print(fig1_handle, 'methods_paper_figures\figure_sc_scatter_1','-dpng','-r300');
print(fig2_handle, 'methods_paper_figures\figure_sc_scatter_2','-dpng','-r300');
print(fig3_handle, 'methods_paper_figures\figure_sc_scatter_3','-dpng','-r300');
print(fig4_handle, 'methods_paper_figures\figure_sc_scatter_4','-dpng','-r300');

%% Suplimental Material Animation.

close all;
clear;
clc;

amino_index = 10;

data_structs = sidechain_atom_extraction();
num_domains = length(data_structs);
input_data = process_sidechain_data(data_structs);


output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);

parfor frame_index = 1:360
    az = frame_index - 1;
    el = 0;

    fig1_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

    ax1 = fig1_handle.CurrentAxes;
    
    axis(ax1, 'off');
    set(ax1, 'DataAspectRatio', [1 1 1]);
    view(ax1,az,el);
    title(ax1,'');
    
    drawnow;
    
    frames(frame_index) = getframe(fig1_handle);
end
fprintf("Done.\n")

videowriter = VideoWriter('methods_paper_figures\sidechain_scatter_rotating_animation.mp4');
videowriter.FrameRate = 10;
open(videowriter);

for i = 1:length(frames)
    frame = frames(i);
    writeVideo(videowriter,frame);
end
close(videowriter);

fprintf("File Written.\n")

delete(gcp('nocreate'));