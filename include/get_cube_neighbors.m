function [neighbors] =get_cube_neighbors(number, cube_size,num_slices)
    neighbors = [];
    
    level = floor(number / (cube_size^2 + 1)) + 1;

    % Add next neighbor in cube row
    if(mod(number, cube_size) > 0)    
        neighbors = union(neighbors, number + 1);
    end
    
    % Add next neighbor in cube columns
    if( number + cube_size <= level*cube_size ^ 2)
        neighbors = union(neighbors, number + cube_size);
    end
    
    % Add next neighbor in cube slices
    if( number + cube_size^2 <= (num_slices * cube_size^2))
        neighbors = union(neighbors, number + cube_size^2);
    end
    
end