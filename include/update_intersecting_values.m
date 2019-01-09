function [] = update_intersecting_values( SubVn_a,SubVn_b, intersection,cube_size)                                                 
    %%
    % Takes two SubVn numbers and merges SubVn_b into SubVn_a if they
    % overlap. 
    % SubVn_a < SubVn_b  (requirement)
    % intersection = 50 (how many pixels intersect with each other)
    % cube_size = 5 (how many subV cubes per side of the whole volume)
    % 
    % Returns a 2D array with the fiber number intersections in 
    % subVn1 and in subVn2
    %%
    %SubVn_a = 7;
    %SubVn_b = 8;
    %intersection=50;
    %cube_size=5;
    
    dir_a = ['SUBVOLUMES/sV' num2str(SubVn_a) '/fibers_info/Vn.mat'];
    dir_b = ['SUBVOLUMES/sV' num2str(SubVn_b) '/fibers_info/Vn.mat'];
   
    
    data_a  = load(dir_a);
    data_b = load(dir_b);
    
    V_a = data_a.Vn;
    V_a_info = data_a.Vn_info;
    
    V_b = data_b.Vn;
    V_b_info = data_b.Vn_info;
    
    [rows,cols,slices] = size(V_a);
    %% Default Starting and Ending Points
    cols_a1 = 1;
    cols_a2 = cols;
    cols_b1 = 1;
    cols_b2 = cols;
        
    rows_a1 = 1;
    rows_a2 = rows;
    rows_b1 = 1;
    rows_b2 = rows;
        
    slcs_a1 = 1;
    slcs_a2 = slices;
    slcs_b1 = 1;
    slcs_b2 = slices;
        
    %% Custom Start and Ending Points
    
    if(SubVn_a == SubVn_b - 1)
        % if b is next column
        cols_a1 = cols - intersection + 1;
        cols_a2 = cols;
        
        cols_b1 = 1;
        cols_b2 = intersection;
    elseif(SubVn_a == SubVn_b - cube_size)
        % if b is next row
        rows_a1 = rows - intersection + 1;
        rows_a2 = rows;
        
        cols_b1 = 1;
        cols_b2 = intersection;
    
    elseif(SubVn_a == SubVn_b - (cube_size^2))
        % if b is next z level
        slcs_a1 = slices - intersection + 1;
        slcs_a2 = slices;
        
        slcs_b1 = 1;
        slcs_b2 = intersection;
    else
        disp('Cubes not touching in straigh line');
    end
    
    
    %% Read/Crop/Multiply Images
        
    V_a_temp = V_a(rows_a1:rows_a2,...
              cols_a1:cols_a2,...
              slcs_a1:slcs_a2);
                   
    V_b_temp = V_b(rows_b1:rows_b2,...
              cols_b1:cols_b2,...
              slcs_b1:slcs_b2);
                
    % Here I do dilatation of the fibers
                
    % Find Fiber Intersection
    SubV_Intersection = V_a_temp .* V_b_temp; 
    
     %% Find Values of Intersection
    overlapping_indices_intersection = find(SubV_Intersection(:) > 0);
    
    values_in_A = unique(V_a_temp(overlapping_indices_intersection));
    values_in_B = unique(V_b_temp(overlapping_indices_intersection));
    
    %values_of_intersection = [values_in_A';...
     %                         values_in_B'];   
    %%
    fibers_merged = [];
    num_fibers_overlaping = length(values_in_A);
    disp(['Found ' num2str(num_fibers_overlaping) ' Fibers Overlapping']);

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
       candidates_unique = unique(candidates_fiber2_raw);
       fiber2 = mode(candidates_fiber2_raw);
       
       angleY1 = V_a_info(fiber1, 4);
       angleZ1 = V_a_info(fiber1, 5);
       
       angleY2 = V_b_info(fiber2, 4);
       angleZ2 = V_b_info(fiber2, 5);
       
       ang_dif_Y = abs(angleY1 - angleY2);
       ang_dif_Z = abs(angleZ1 - angleZ2);
       
       while(~ismember(fiber2, values_in_B) && fiber2 > 0)
           candidates_fiber2_raw(find(candidates_fiber2_raw == fiber2)) = [];
           fiber2 = mode(candidates_fiber2_raw);
           
           angleY2 = V_b_info(fiber2, 4);
           angleZ2 = V_b_info(fiber2, 5);
           
           ang_dif_Y = abs(angleY1 - angleY2);
           ang_dif_Z = abs(angleZ1 - angleZ2);
       
       end
       if(fiber2 == 0 || isnan(fiber2))
           continue;
       end
       
       % Remove Fiber from Possible Connections in Volume B
       values_in_B(find(values_in_B == fiber2)) = [];
       
       % Update Subvolume 2
       indx_to_update = find(V_b == fiber2);
       V_b(indx_to_update) = fiber1;
       
       %TO DO: Update text file
       % Merge Fiber 1 and Fiber 2 id:
       V_b_info(fiber2, 1) = V_a_info(fiber1, 1);
       
       % Merge raddii
       new_R = (V_a_info(fiber1, 2) + V_b_info(fiber2, 2))/2;
       V_a_info(fiber1, 2) = new_R;
       V_b_info(fiber2, 2) = new_R;
       
       % Merge length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TWEAK
       new_H = V_a_info(fiber1, 3) + V_b_info(fiber2, 3);
       V_a_info(fiber1, 3) = new_H;
       V_b_info(fiber2, 3) = new_H;
       
       % Merge x,y,z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TWEAK
       new_X = (V_a_info(fiber1, 6) + V_b_info(fiber2, 6))/2;
       new_Y = (V_a_info(fiber1, 7) + V_b_info(fiber2, 7))/2;
       new_Z = (V_a_info(fiber1, 8) + V_b_info(fiber2, 8))/2;
       V_a_info(fiber1, 6) = new_X;
       V_b_info(fiber2, 6) = new_X;
       
       V_a_info(fiber1, 7) = new_Y;
       V_b_info(fiber2, 7) = new_Y;
       
       V_a_info(fiber1, 8) = new_Z;
       V_b_info(fiber2, 8) = new_Z;
       
       
       disp(['Merged Fiber ' num2str(fiber1) ' and ' num2str(fiber2)]);
       
       
       %debug_volume(fiber1_Volume, V_b_temp,fiber1_intersections,'DEBUG');
    end
    
    % Save Updated Results
    Vn = V_b;
    Vn_info = V_b_info;
    save(dir_b,'Vn', 'Vn_info');
end