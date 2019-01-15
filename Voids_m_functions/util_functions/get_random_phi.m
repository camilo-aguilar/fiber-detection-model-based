function [ phi,expansion_weight] = get_random_phi(r_gdc,mpp_i_n,dims)

    % r_gdc1 :1 circle, 2 square ,3 circle array ,4 chess board
    % r_gdc2 :-1 expand, 1 contract
    % r_gdc3 : radious 
    % r_gdc4 : center for shape type 1
    % r_gdc5 : distance for each radious type 3 and 4
    
    
    if(mpp_i_n == 1)
        %Large Rectangle to contract
        phi = ac_SDF_2D('rectangle', dims, 5);
        expansion_weight = -1;
    elseif(mpp_i_n == 2)
        %Small Rectangle to expand
        phi = ac_SDF_2D('rectangle', dims, (min(dims(:)))/2 - 20);
        expansion_weight = 1;
    elseif(mpp_i_n == 3)
        %Circles to contract
        phi = ac_SDF_2D('circle_array', dims, 28 , 60);
        expansion_weight = -1;
    elseif(mpp_i_n == 4)
        %Circles to expand
        phi = ac_SDF_2D('circle_array', dims, 20 , 100);
        expansion_weight = 1;
    elseif(mpp_i_n == -1)
        rd = 890;
        center(1) = round(dims(1)/2)+10;
        center(2) = round(dims(2)/2)-60;
        expansion_weight = -1;
        phi = ac_SDF_2D('circle', dims, center, rd);  
        phi(phi<0) = -inf;
        
    else
        %Circle to contract
        rd = round(r_gdc(1)*90)+30;
        center(1) = round( r_gdc(2)*(dims(1)-rd)+ rd);
        center(2) = round( r_gdc(3)*(dims(2)-rd) + rd);
        expansion_weight = -1;
        phi = ac_SDF_2D('circle', dims, center, rd);  
        phi(phi<0) = -inf;
    end

    return;
end

