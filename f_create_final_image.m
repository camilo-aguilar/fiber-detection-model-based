%% File: Create Final Image
%  Author: Camilo Aguilar
%  Input:  Directory 'SUBVOLUMES/sVn' (where n = 1,2,3,...75)
%          Each folder 'SUBVOLUMES/sVn' has 5 subfolders to store
%          information: data,fibers,seg,voids,fibers_info
%
%  Output: FINAL_RESULT/final_n.tiff
%          final_n.tiff is an rgb image
%          where fiber_number = r*256^2 + g*256 + b
%%

tic
batch_size = 450;
cubes_per_batch = 5;
windows_per_level = 25;
num_levels = 3;
pixels_overlap = 50;

coordinates_start = 1:400:1601;
coordinates_end = 450:400:2050;


directory = 'FINAL_RESULT';

% % Create Images
if(~exist(directory,'dir'))
    mkdir(directory);
end

directory = [directory '/'];

for i = 1:450*3
     if(i < 10)
         name = [directory 'final_000' num2str(i) '.tif'];
     elseif(i < 100)
         name = [directory 'final_00' num2str(i) '.tif'];

     elseif(i < 1000)
         name = [directory 'final_0' num2str(i) '.tif'];
     else
         name = [directory 'final_' num2str(i) '.tif'];
     end 
     im_rgb = zeros(2050,2050,3, 'uint8');
     imwrite(im_rgb, name);
end

disp('Getting Coordinates');
coordinates(1) = struct('rows', 1:batch_size, 'cols', 1:batch_size, 'slices', 1:batch_size);
cube_number = 1;
for zn=1:num_levels
    for in=1:sqrt(windows_per_level)
        for jn=1:sqrt(windows_per_level)          
            x = [coordinates_start(in), coordinates_end(in)];
            y = [coordinates_start(jn),coordinates_end(jn)];
            z = [coordinates_start(zn),coordinates_end(zn)];
            coordinates(cube_number) = struct('rows', x, 'cols', y, 'slices', z);
            cube_number = cube_number + 1;
        end
            
    end
end


disp('Updating Volume');

fileID = fopen('FINAL_RESULTS.txt','w');

for cube_number=1:75
    cube_coordinates = coordinates(cube_number);
    dir_a = ['SUBVOLUMES/sV' num2str(cube_number) '/fibers_info/Vn.mat'];
    data_a = load(dir_a);
    V_a_info = data_a.Vn_info;
    V_a = data_a.Vn;
    fprintf(fileID,'%4.4f,%4.4f,%4.4f,%4.4f,%4.4f,%4.4f,%4.4f,%4.4f\n',V_a_info);

    for slc=cube_coordinates.slices(1):cube_coordinates.slices(2)
        counter = slc - cube_coordinates.slices(1) + 1;
        
        if(slc < 9)
          name = [directory 'final_000' num2str(slc) '.tif'];
        elseif(slc < 99)
          name = [directory 'final_00' num2str(slc) '.tif'];
        elseif(slc < 999)
          name = [directory 'final_0' num2str(slc) '.tif'];
        else
          name = [directory 'final_' num2str(slc) '.tif'];
        end
        
        % Get specific coordinates of interest
        i_coord = cube_coordinates.rows(1):cube_coordinates.rows(2); 
        j_coord = cube_coordinates.cols(1):cube_coordinates.cols(2);
        
        % Read image and crop it
        im_rgb = imread(name);
        im_rgb_d = double(im_rgb(i_coord, j_coord,:));
        
        old_values = im_rgb_d(:,:,1)*256^2 + im_rgb_d(:,:,2)*256 + im_rgb_d(:,:,3);
        new_values = V_a(:,:,counter);
        
        % Merge Overlapping Regions together
        updated_values = max(new_values,old_values);
        
        % Calculate red component (third bit)
        r = uint8(bitand(updated_values, bitshift(255, 16)));
        
        % Calculate green component (second bit)
        g = uint8(bitand(updated_values, bitshift(255, 8)));
        
        % Calculate blue component (first bit)
        b = uint8(bitand(updated_values, 255));
        
        
        % Update Image
       im_rgb(i_coord, j_coord,1) = r;
       
       im_rgb(i_coord, j_coord,2) = g;
       
       im_rgb(i_coord, j_coord,3) = b;
       
       % Save Image
       imwrite(im_rgb,name);
       
    end
    
    disp(['Updating Volume' num2str(cube_number)]);
end
fclose(fileID);

toc
