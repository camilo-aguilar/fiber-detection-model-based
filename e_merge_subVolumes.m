%% File: merge_subvolumes
%  Author: Camilo Aguilar
%  Function: Merge Overlapping Fibers
%
%  Input:  Files 'SUBVOLUMES/sV#/fibers_info/Vn.mat' (where # = 1,2,3,...75)
%  Output: 
%          'SUBVOLUMES/sV#/fibers_info/Vn.mat' but with overlapping fibers
%           merged
addpath('include');

%% Merge Pre-Procced-subVolumes pairwise
cube_size = 5;
blocks_per_height = 3;
pixel_overlap = 50;
ANGLE_THRESHOLD = 12;

 tic
 merged_fibers = 0;

 for SubVn_a=1:75
     disp(['Merging Volume: ' num2str(SubVn_a)]);
     % Load Volume Vo to memory
     dir_a = ['SUBVOLUMES/sV' num2str(SubVn_a) '/fibers_info/Vn.mat'];
     data_a  = load(dir_a);
     V_a = data_a.Vn;
     V_a_info = data_a.Vn_info;
     offset_a = data_a.offset;
     
     % Get dimensions
     [rows,cols,slices] = size(V_a);
    
     % Get all the adjacent neighbors to Vo (in raster order)
     neighbors_list = get_cube_neighbors(SubVn_a, cube_size, blocks_per_height);
     
     %% UPDATE ADJACENT VOLUMES
     for Vf_n=neighbors_list
        % Load Volume B
        dir_b = ['SUBVOLUMES/sV' num2str(Vf_n) '/fibers_info/Vn.mat'];
        data_b = load(dir_b);
        V_b = data_b.Vn;
        V_b_info = data_b.Vn_info;
        offset_b = data_b.offset;
    
        [rows_a, cols_a, slices_a, rows_b, cols_b, slices_b] = get_indexes(SubVn_a, Vf_n,rows,cols,slices,pixel_overlap,cube_size);
        
        
        %% Crop/Multiply Images
        V_a_temp = V_a(rows_a, cols_a, slices_a);
        V_b_temp = V_b(rows_b, cols_b, slices_b);
        
        
        % Here I could dilatation of fibers
                
        % Find Fiber Intersection
        SubV_Intersection = V_a_temp .* V_b_temp; 
        
        %% Find Values of Intersection
        overlapping_indices_intersection = find(SubV_Intersection(:) > 0);
    
        values_in_A = unique(V_a_temp(overlapping_indices_intersection));
        values_in_B = unique(V_b_temp(overlapping_indices_intersection));
          
        %%%%%%%%%%%%%% update_intersecting_values( Vo_n,Vf_n,pixel_overlap,blocks_per_row); 
        
        num_fibers_overlaping = length(values_in_A);
        disp(['Found ' num2str(num_fibers_overlaping) ' Fibers Overlapping']);
        merged_counter = 0;
        merged_values_b = [];
        for indx=1:num_fibers_overlaping
            fiber1 = values_in_A(indx);
       
            % Find all values in V1 temp where this fiber is located
            fiber1_indx = find(V_a_temp == fiber1);
       
            % Create a empty volume only with this Fiber inside
            fiber1_Volume = zeros(size(V_a_temp));
            fiber1_Volume(fiber1_indx) = 1;
       
            % Filter SubVn2 with intersections
            fiber1_intersections = fiber1_Volume .* V_b_temp;
            intersection2 = find(fiber1_intersections > 0);
       
            % Find best fiber to merge
            candidates_fiber2_raw = V_b_temp(intersection2);
       
            % Check Candidates Angles Coherence and Area of Intersection
            fiber2 = mode(candidates_fiber2_raw);

            %Check if Fiber should be merged with two Fibers inside the same
            %volume
            fiber2_offset = find(V_b_info(1,:) == fiber2);
            if(length(fiber2_offset) > 1)
                V_b_info(2:end,fiber2_offset(1)) = mean(V_b_info(2:end,fiber2_offset),2);
                V_b_info(:,fiber2_offset(2:end)) = [];
                fiber2_offset = fiber2_offset(1);
                disp('Inner Merge in V2');
            end
            
            fiber1_offset = find(V_a_info(1 ,:) == fiber1);
            if(length(fiber1_offset) > 1)
                V_a_info(2:end,fiber1_offset(1)) = mean(V_a_info(2:end,fiber1_offset),2);
                V_a_info(:,fiber1_offset(2:end)) = [];
                fiber1_offset = fiber1_offset(1);
                disp('Inner Merge in V1');
            end
            
            
            angleY1 = V_a_info(4,fiber1_offset);
            angleZ1 = V_a_info(5,fiber1_offset);
       
            angleY2 = V_b_info(4,fiber2_offset);
            angleZ2 = V_b_info(5, fiber2_offset);
            
            if(abs(360 - angleY1) < ANGLE_THRESHOLD)
                angleY1 = abs(360 - angleY1);
            end
            
            if(abs(360 - angleY2) < ANGLE_THRESHOLD)
                angleY2 = abs(360 - angleY2);
            end
            
            if(abs(360 - angleZ1) < ANGLE_THRESHOLD)
                angleZ1 = abs(360 - angleZ1);
            end

            if(abs(360 - angleZ2) < ANGLE_THRESHOLD)
                angleZ2 = abs(360 - angleZ2);
            end

            ang_dif_Y = abs(angleY1 - angleY2);
            ang_dif_Z = abs(angleZ1 - angleZ2);
            
            if((ismember(fiber2, merged_values_b)))
                ang_dif_Y =100;
                ang_dif_Z = 100;
            end    
            % HERE IS THE UNION. IF 
            while((ismember(fiber2, merged_values_b)...
                   || ang_dif_Y > ANGLE_THRESHOLD ...
                   || ang_dif_Z > ANGLE_THRESHOLD) ...
                   && fiber2 > 0)
                candidates_fiber2_raw(find(candidates_fiber2_raw == fiber2)) = [];
                fiber2 = mode(candidates_fiber2_raw);
                
                if(fiber2 > 0)
                    fiber2_offset = find(V_b_info(1,:) == fiber2);
                    angleY2 = V_b_info(4, fiber2_offset);
                    angleZ2 = V_b_info(5, fiber2_offset);
                    
                    if(abs(360 - angleY2) < ANGLE_THRESHOLD)
                        angleY2 = abs(360 - angleY2);
                    end
                     
                     if(abs(360 - angleZ2) < ANGLE_THRESHOLD)
                        angleZ2 = abs(360 - angleZ2);
                     end
                    
                    
                    ang_dif_Y = abs(angleY1 - angleY2);
                    ang_dif_Z = abs(angleZ1 - angleZ2);
                else
                    break;
                end
           
                
       
            end
            if(fiber2 == 0 || isnan(fiber2))
                continue;
            end
       

            % Remove Fiber from Possible Connections in Volume B
            merged_values_b = union(merged_values_b, fiber2);
       
            % Update Subvolumes
            % find lowest fiber:
            [~, largest] = max([fiber1, fiber2]);
            
            if(largest == 1)
                % Replace Fiber # in volume 1
                indx_to_update = find(V_a == fiber1);
                V_a(indx_to_update) = fiber2;
                V_a_info(1 ,fiber1_offset) = V_b_info(1, fiber2_offset);
                
            else
                % Replace Fiber # in volume 2
               indx_to_update = find(V_b == fiber2);
               V_b(indx_to_update) = fiber1;
               V_b_info(1 ,fiber2_offset) = V_a_info(1, fiber1_offset);
            end
            
           

           % Merge raddii
           new_R = (V_a_info(2, fiber1_offset) + V_b_info(2,fiber2_offset))/2;
           V_a_info(2, fiber1_offset) = new_R;
           V_b_info(2, fiber2_offset) = new_R;

           % Merge length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TWEAK
           new_H = V_a_info(3, fiber1_offset) + V_b_info(3, fiber2_offset);
           V_a_info(3, fiber1_offset) = new_H;
           V_b_info(3, fiber2_offset) = new_H;

           % Merge x,y,z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TWEAK
           new_X = (V_a_info(6, fiber1_offset) + V_b_info(6,fiber2_offset))/2;
           new_Y = (V_a_info(7, fiber1_offset) + V_b_info(7,fiber2_offset))/2;
           new_Z = (V_a_info(8, fiber1_offset) + V_b_info(8,fiber2_offset))/2;
           V_a_info(6,fiber1_offset) = new_X;
           V_b_info(6,fiber2_offset) = new_X;

           V_a_info(7,fiber1_offset) = new_Y;
           V_b_info(7,fiber2_offset) = new_Y;

           V_a_info(8,fiber1_offset) = new_Z;
           V_b_info(8,fiber2_offset) = new_Z;

           merged_counter = merged_counter + 1;
           disp(['Merged Fibers: ' num2str(fiber1) ' and ' num2str(fiber2)]);
           merged_fibers = merged_fibers + 1;
        end
        disp(['Done with neighbor ' num2str(Vf_n) ',merged: ' num2str(merged_counter) ]);
        % Save Updated Results
        Vn = V_a;
        Vn_info = V_a_info;
        offset = offset_a;
        save(dir_a,'Vn', 'Vn_info', 'offset');
        
        Vn = V_b;
        Vn_info = V_b_info;
        offset = offset_b;
        save(dir_b,'Vn', 'Vn_info', 'offset');
        
     end
 end
toc
