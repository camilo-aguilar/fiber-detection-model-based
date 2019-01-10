%% File: Crop Original Volume
%  Author: Camilo Aguilar
%  Function: Crop the original volume into 2050 x 2050 x 1350 volume
%
%  Input:  Original uint16 images located in folder "input_dir"
%  Output: Cropped uint8 images located in folder "CROPPED_FILES"
%
%%
%% Set parameters
tic
input_dir = 'INPUT_FILES/';
output_dir = 'CROPPED_FILES';

if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end
num_images = 1392;
output_dir = [output_dir '/'];

%% Start of the script
parfor i=1:num_images
    if i > 999
        number = num2str(i);
    elseif i > 99
        number = ['0' num2str(i)];
    elseif i > 9
        number = ['00' num2str(i)];
    else
        number = ['000' num2str(i)];
    end
    
    name = ['STACK' number];
    im = imread([input_dir name '.tif']);
    
    im = im(200:end-311,280:end-231);
    imc = uint8(round(double(im) / 2^8));
    imwrite(imc, [output_dir name '_cropped.tif']);
    
end
toc 