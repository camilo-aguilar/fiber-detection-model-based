original = zeros(450,450,3,450);

for i=1:450
original(:,:,:,i) = imread(['mer_' num2str(i) '.tif']);
end

transformed = zeros(450,450,450);
for i=1:450
    transformed(:,:,i) = original(:,:,3, i) + (256.0 * original(:,:,2,i)) + (256.0^2 * original(:,:,1,i));
end

fileID = fopen('merged_fiber.txt','r');
    Vn_info = fscanf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f',[8 Inf]);
fclose(fileID);

single = zeros(size(transformed));
single(find(transformed == 68)) = 1;
debug_volume( single,single,single, 'DEBUG' )