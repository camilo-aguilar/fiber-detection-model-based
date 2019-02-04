%% File: script1
%  Author: Camilo Aguilar
%%

% Crop the original volume into 2050 x 2050 x 1350 volume
% a_crop_original_volume

% Split the 2050 x 2050 x 1350 volume into 75 overlapping cubes
b_split_subVolumes

% Send 75 scripts to cluster to find fibers using C program
c_find_fibers_voids

