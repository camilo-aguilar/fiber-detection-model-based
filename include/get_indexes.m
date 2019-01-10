function [rows_a,cols_a, slices_a, rows_b, cols_b,slices_b] = get_indexes( SubVn_a, SubVn_b, rows,cols,slices,intersection,cube_size )
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
        
        rows_b1 = 1;
        rows_b2 = intersection;
    
    elseif(SubVn_a == SubVn_b - (cube_size^2))
        % if b is next z level
        slcs_a1 = slices - intersection + 1;
        slcs_a2 = slices;
        
        slcs_b1 = 1;
        slcs_b2 = intersection;
    else
        disp('Cubes not touching in straigh line');
    end
    
    rows_a = rows_a1:rows_a2;
    cols_a = cols_a1:cols_a2;
    slices_a = slcs_a1:slcs_a2;
    
    rows_b = rows_b1:rows_b2;
    cols_b = cols_b1:cols_b2;
    slices_b = slcs_b1:slcs_b2;
    

end

