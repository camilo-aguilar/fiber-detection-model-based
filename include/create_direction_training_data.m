fibers_rgb = zeros(450,450,3,450);
voids = zeros(450,450,450,'uint8');

for i=1:450
fibers_rgb(:,:,:,i) = imread(['mer_' num2str(i) '.tif']);
voids(:,:,:,i) = imread(['../filtered_voids/filtered_void' num2str(i) '.tif']);
end

transformed = zeros(450,450,450);
for i=1:450
    transformed(:,:,i) = fibers_rgb(:,:,3, i) + (256.0 * fibers_rgb(:,:,2,i)) + (256.0^2 * fibers_rgb(:,:,1,i));
end

txt_dir = '../fibers_info/merged_fiber.txt';
V = get_subVolume( img_dir, 450);
fileID = fopen(txt_dir,'r');
Vn_info = fscanf(fileID,'%d,%f,%f,%f,%f,%f,%f,%f',[8 Inf]);
Vn_info = Vn_info(:,2:end);
Vn_info(1,:) = Vn_info(1,:) - 1;
fclose(fileID);

single = zeros(size(transformed), 'uint8');

fiber_classes = zeros(3,3,'uint8');
counter = 0;
for ay=1:3
    for az=1:3
      fiber_classes(ay,az) = counter;
      counter = counter + 1;
    end
end

for k=1:450
    for i=1:450
        for j=1:450
            fiber_n = transformed(i,j,k);
            if(voids(i,j,k) > 0)
                single(i,j,k) = 250;
            end
            if(fiber_n > 0)
                angleY = Vn_info(fiber_n,4) * 2;
                angleZ = Vn_info(fiber_n,5) * 2;
                
                indy = floor(angleY/60.001) + 1;
                indz = floor(angle/60.001) + 1;
                single(i,j,k) = 25 * fiber_classes(indy,indz);
            end
        end
    end
end

for i=1:450
    imwrite(single(:,:,i),['debug_' num2str(i) '.jpeg']);
end
