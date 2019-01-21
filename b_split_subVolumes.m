%% File: split_subvolumes
%  Author: Camilo Aguilar
%  Function: Split the 2050 x 2050 x 1350 volume into 75 overlapping cubes
%
%  Input:  Original images located in folder "INPUT_FILES/"
%  Output: Directory 'SUBVOLUMES/sVn' (where n = 1,2,3,...75)
%          Each folder 'SUBVOLUMES/sVn' has 5 subfolders to store
%          information: data,fibers,seg,voids,fibers_info
%          This sccript only populates data and seg

addpath('include');

if(~exist('SUBVOLUMES', 'dir' ))
    mkdir('SUBVOLUMES')
end
%% Set parameters
input_dir = 'INPUT_FILES/';

% These parameters are fixed by now
VOLUME_SIZE = [2050, 2050, 1392];
subvolume_size = 450;
subvolumes_per_level = 25;
num_levels = 3;



%% Create Splitting Coordiates Table
% Creates 75 BLOCKS OF (5 x 5 x 3) SUBVOLUMES EACH CONTAINING 450^3 voxels 
% Each bloc overlaps (50x475x475) voxels with a neighbot in each side
% Divides in raster order: 
%   left to right (columnwise)
%   top to bottom (rowise)
%   top to bottom (slices)
coordinates_start = 1:400:1601;
coordinates_end = 450:400:2050;

coordinates(1) = struct('rows', 1:subvolume_size, 'cols', 1:subvolume_size, 'slices', 1:subvolume_size);
cube_number = 1;
for zn=1:num_levels
    for in=1:sqrt(subvolumes_per_level)
        for jn=1:sqrt(subvolumes_per_level)
            
            x = (coordinates_start(in):coordinates_end(in)) + 200;
            y = (coordinates_start(jn):coordinates_end(jn)) + 280;
            z = coordinates_start(zn):coordinates_end(zn);
            coordinates(cube_number) = struct('rows', x, 'cols', y, 'slices', z);
            cube_number = cube_number + 1;
        end          
    end
end

%% Create directories in each SubVolume
% Data:         Cropped Subvolume
% fibers:       Detected Fibers
% Voids:        Detected Voids
% seg:          EM/MPM segmentation
% fibers_info:  Text files and .mat files with processed information

disp('Starting to Process SubVolumes');
parfor cube_number=1:75
    disp(['Starting Volume: ' num2str(cube_number)]);
    disp('Creating Diretories');
    dir_name = ['SUBVOLUMES/sV' num2str(cube_number)];
    % Create Folders
    mkdir(dir_name)
    mkdir([dir_name '/data'])
    mkdir([dir_name '/fibers'])
    mkdir([dir_name '/voids'])
    mkdir([dir_name '/seg'])
    mkdir([dir_name '/fibers_info'])
    
    disp('Creating Cropped SubVolume...');
    % Save Subvolume
    locations = coordinates(cube_number);
    V = crop_subVolume(locations.rows,locations.cols,locations.slices,input_dir);
    
    disp('Saving and Segmenting Cropped SubVolume...');
    for i=1:size(V,3)
       dir_to_write = [dir_name '/data/'];
       dir_to_write_seg = [dir_name '/seg/'];
       if(i > 999)
           name_original = [dir_to_write 'subV_' num2str(i) '.tif'];
           name_segmented = [dir_to_write_seg 'subV_seg_' num2str(i) '.tif'];
       elseif(i > 99)
           name_original = [dir_to_write 'subV_0' num2str(i) '.tif'];
           name_segmented = [dir_to_write_seg 'subV_seg_0' num2str(i) '.tif'];
       elseif(i > 9)
           name_original = [dir_to_write 'subV_00' num2str(i) '.tif'];
           name_segmented = [dir_to_write_seg 'subV_seg_00' num2str(i) '.tif'];
       else
           name_original = [dir_to_write 'subV_000' num2str(i) '.tif'];
           name_segmented = [dir_to_write_seg 'subV_seg_000' num2str(i) '.tif'];
       end
       segged = emmpm(uint8(V(:,:,i)),6);
       segged = uint8(round(255 * double(segged)/6));
       imwrite(uint8(V(:,:,i)), name_original);
       imwrite(segged, name_segmented);
       
    end
    
end

%%
