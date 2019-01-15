function [Ext] = Snake_Energy(rx,ry,energy_image,g)
    %% External Energy  
    if(isempty(rx))
        Ext = inf;
        return
    end
    
    segged_I = DrawSegmentedArea2D(zeros(size(energy_image)),rx,ry);
    segged_I = segged_I';
    
    Inner_black_image = (segged_I.* energy_image);
    Inner_Area = length(find(segged_I > 0));
    
    
    
    [~,false_black_area_i] = find(Inner_black_image(:) < 0);
    Inner_black_image(false_black_area_i) = Inner_black_image(false_black_area_i)*10;
    
    Inner_black_area = sum(Inner_black_image (:));
    
    %if(E_black ==0)
    %    Ext = inf;
    %    return;
    %end
    
    E_edge = 0;
    if(exist('g','var'))
        edge_values = interp2(g,rx,ry);
        [~,false_edge_i] = find(edge_values <= 0);
        edge_values(false_edge_i) = edge_values(false_edge_i) - 5; 
        E_edge = sum(edge_values);
        
    end
    
    Ext = Inner_black_area + 4*E_edge./length(rx);
    if(Inner_black_area <= 0 || E_edge < 0)
        Ext = inf;
    end


    %%
end