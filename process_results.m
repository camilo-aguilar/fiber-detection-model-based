%% File: Process Results
%  Author:  Camilo Aguilar
%  Function: Read final results from void and fiber detection
%
%  Input:  RGB Images 'FINAL_RESULTS/final_i.tiff' (i = 1:1350)
%          Binary Images 'SEGMENTED_VOIDS/voids_n.tiff' (i = 1:1350)
%          Text File: FINAL_RESULTS.txt
%          
%  Output: 
%          
%
%%

fiber_pixels = 0;
void_pixels = 0;
matrix_pixels =0;

% Dimensions of Input/Output Image
rows = 2050;
cols = 2050;

% To count only items inside the circular sample
[x,y] = meshgrid(-rows/2+1:rows/2,-cols/2+1:cols/2);
d = 910.^2 - (x+10).^2 - (y+5).^2 ; 
circle = sign(d);
sample_pixels = length(find(circle > 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% In case that want to store final volume (DEBUG)
%final_volume = zeros(rows,cols,150,'uint16');
%final_v_counter = 1;
%vol_num = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 for num=1:1350
     disp(num);
     if num < 10
         fibers = imread(['FINAL_RESULT/final_000' num2str(num) '.tif']);
         voids = imread(['SEGMENTED_VOIDS/voids_000' num2str(num) '.tiff']);
     elseif num < 100
         fibers = imread(['FINAL_RESULT/final_00' num2str(num) '.tif']);
         voids = imread(['SEGMENTED_VOIDS/voids_00' num2str(num) '.tiff']);
     elseif num < 1000
         fibers = imread(['FINAL_RESULT/final_0' num2str(num) '.tif']);
         voids = imread(['SEGMENTED_VOIDS/voids_0' num2str(num) '.tiff']);
     else
         fibers = imread(['FINAL_RESULT/final_' num2str(num) '.tif']);
         voids = imread(['SEGMENTED_VOIDS/voids_' num2str(num) '.tiff']);
     end
     
     [rows,cols,channels] = size(fibers);
     
     % Transform fibers to uint16
     fibers = double(fibers);
     fibers_non_zero = fibers(:,:,3) + (256.0 * fibers(:,:,2)) + (256.0^2 * fibers(:,:,1));
     fibers_non_zero = uint16(fibers_non_zero);
     
     % Transform voids to uint16
     voids = (uint16(sign(voids))*65535);
     
     % Merge Fibers and Voids by setting voids in overlapping indices to 0
     overlap_indices = find((fibers_non_zero .* voids) > 0);
     voids(overlap_indices) = 0;
     
     % Find number of fibers/voids in current slice
     current_fiber_pixels = length(find(fibers_non_zero > 0));
     current_void_pixels = length(find(voids > 0));
     
     % Update counters for matrix,fiber and void pixels
     fiber_pixels = fiber_pixels + current_fiber_pixels;
     void_pixels = void_pixels + current_void_pixels; 
     matrix_pixels =  matrix_pixels + (sample_pixels  - current_fiber_pixels - current_void_pixels);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %final_volume(:,:,final_v_counter) = fibers_non_zero + voids;
     %final_v_counter = final_v_counter + 1;
     %if(final_v_counter == 150)
     %    final_v_counter = 1;
     %    save(['finalV' num2str(vol_num)],'final_volume');
     %    disp('Saved');
     %    disp(fiber_pixels/matrix_pixels);
     %    vol_num = vol_num + 1;
     %end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end

%% Get Statistics
fileID = fopen('FINAL_RESULTS.txt','r');
Vn_info = fscanf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f',[8 Inf]);
fclose(fileID);

% Process Statistics
num_entries = length(Vn_info);

Fiber_lengths = zeros(1,num_entries);
Fiber_radii = zeros(1,num_entries);
Angle_z_axis = zeros(1,num_entries);
Angle_y_axis = zeros(1,num_entries);

%Fill out statistics (Merged Fibers preserve the averages)
for i=1:num_entries
   Fiber_radii(Vn_info(1,i)) = Vn_info(2,i);
   Fiber_lengths(Vn_info(1,i)) = Vn_info(3,i);   
   Angle_z_axis(Vn_info(1,i)) = Vn_info(4,i);
   Angle_y_axis(Vn_info(1,i)) = Vn_info(5,i);
end

threshold = 1;
% Filter Fibers too short/thin
Filtered_Lengths = Fiber_lengths(Fiber_lengths > threshold * 20);
Filtered_Radii = Fiber_radii(Fiber_lengths > threshold);
Filtered_Angle_y = Angle_y_axis(Fiber_lengths > threshold);
Filtered_Angle_z = Angle_z_axis(Fiber_lengths > threshold);


% Create Histogram
% histogram(Filtered_Lengths);

%% SAVE ALL RESULTS
save('Processed_Merged_Results',...
     'fiber_pixels', ...
      'void_pixels', ...
      'matrix_pixels', ...
      'Filtered_Lengths', ...
      'Filtered_Radii', ...
      'Filtered_Angle_y', ...
      'Filtered_Angle_z');
subplot(2,2,1);
histogram(Filtered_Lengths); title('Lengths');
subplot(2,2,2);
histogram(Filtered_Radii); title('Radii');
subplot(2,2,3);
histogram(Filtered_Angle_y); title('Orientation y');
subplot(2,2,4);
histogram(Filtered_Angle_z); title('Orientation z');

fiber_matrix_r = fiber_pixels/matrix_pixels

void_matrix_r = void_pixels/matrix_pixels

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL TO SAVE/SEE COMBINED RESULTS

% Create rgb arrays of 1000 random colors
r_a = floor(255 * rand(1,1000) + 1);
g_a = floor(255 * rand(1,1000) + 1);
b_a = floor(255 * rand(1,1000) + 1);

% Set rgb(1) = 0 for background
r_a(1) = 0;
g_a(1) = 0;
b_a(1) = 0;

% Select which image to save (it is set to random but can choose range of #s)
num = round(600*rand(1,1) + 1)+500;

if num < 10
    fibers = imread(['FINAL_RESULT/final_000' num2str(num) '.tif']);
    voids = imread(['CONTOUR_RESULTS/voids_casselle_' num2str(num) '.tiff']);
elseif num < 100
    fibers = imread(['FINAL_RESULT/final_00' num2str(num) '.tif']);
    voids = imread(['CONTOUR_RESULTS/voids_casselle_' num2str(num) '.tiff']);
elseif num < 1000
    fibers = imread(['FINAL_RESULT/final_0' num2str(num) '.tif']);
    voids = imread(['CONTOUR_RESULTS/voids_casselle_' num2str(num) '.tiff']);
else
    fibers = imread(['FINAL_RESULT/final_' num2str(num) '.tif']);
    voids = imread(['CONTOUR_RESULTS/voids_casselle_' num2str(num) '.tiff']);
end

% Find Fiber Number
r_fi = double(fibers(:,:,1)) * 256^2;
g_fi = double(fibers(:,:,2)) * 256;p
b_fi = double(fibers(:,:,3));
fiber_num = r_fi+g_fi+b_fi;

% Remap Fibers to random colors
index = mod(fiber_num,1000) + 1;
r = uint8(r_a(index));
g = uint8(g_a(index));
b = uint8(b_a(index));

% Filter Fibers to only see the ones inside the sample (filter artifacts)
sample = uint8(max(circle,0));
fiber_color = cat(3,r .* sample,g .* sample, b.* sample);

% Merge with void results (voids = red)
fiber_color(:,:,1) = max(fiber_color(:,:,1),voids(201:2250, 281:2330));

% Show/Save image
imshow(fiber_color);
