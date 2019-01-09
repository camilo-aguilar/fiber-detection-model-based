function [  ] = debug_volume( V1,V2,V3, Vdir )
[rows,cols,slcs] = size(V1);

cV = zeros(rows,cols,3);
for i=1:slcs
    temp1 = V1(:,:,i);
    temp1 = uint8(255 * sign(temp1));
    cV(:,:,1) = temp1;
    
    temp2 = V2(:,:,i);
    temp2 = uint8(255 * sign(temp2));
    cV(:,:,3) = temp2;
    
    temp3 = V3(:,:,i);
    temp3 = uint8(255 * sign(temp3));
    cV(:,:,2) = temp3;
    
    imwrite(cV,[Vdir '/debug_' num2str(i) '.jpeg']);
end

