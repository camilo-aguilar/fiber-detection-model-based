%% File: script2
%  Author: Camilo Aguilar
%%

% Assign a number to each detected fiber and store it
d_preprocess_subVolumes

% Merge Overlapping Fibers
e_merge_subVolumes

% Create Final Image and text file with Merged Fibers 
f_create_final_image