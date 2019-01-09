function [ V ] = get_subVolume( img_dir, num_ims )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    % Directory
    im = imread([img_dir num2str(1) '.tif']);
    [rows,cols, channels] = size(im);
    
    V = zeros(rows,cols,channels, num_ims);
    
    for i=1:num_ims
        V(:,:,:,i) = imread([img_dir num2str(i) '.tif']);
    end
    
    V = squeeze(V);
end

