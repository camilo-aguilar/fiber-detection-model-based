function [ V ] = crop_subVolume( rows_v, cols_v, slices_v, input_dir)
    %%
    % INPUT: 
    %       rows_v,cols_v,slices_v: Coordinate vectors to crop
    % OUTPUT:
    %       V: Sub Volume of Cropped STACK Images at input vectors
    % 
    %%
    if(~exist(input_dir, 'dir'))
        input_dir = 'CROPPED_FILES/';
    end
    V = zeros(length(rows_v), length(cols_v), length(slices_v));
    counter = 1;
    for i=slices_v
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
        im = im(rows_v,cols_v);
        im = im / 2 ^ 8;
        V(:,:,counter) = im;
        counter = counter + 1;
    end


end

