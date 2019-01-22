function [ ] = voids_2D( Num_z)
%
%
%
    addpath('Voids_m_functions/util_functions');

    tic


    show_result = 0;

    %% Input Image
    if(Num_z) >999
        I = imread(['INPUT_FILES/STACK' num2str(Num_z) '.tif']);
    elseif(Num_z) >99
        I = imread(['INPUT_FILES/STACK0' num2str(Num_z) '.tif']);
    elseif(Num_z) >9
        I = imread(['INPUT_FILES/STACK00' num2str(Num_z) '.tif']);
    else
        I = imread(['INPUT_FILES/STACK000' num2str(Num_z) '.tif']);
    end
    I = double(I)./65535 * 255;

    %[rows,cols] = size(I);

    %[x,y] = meshgrid(-rows/2+1:rows/2,-cols/2+1:cols/2);
    %d = 950.^2 - (x-20).^2 - (y+65).^2 ; 

    %mask = zeros(rows,cols);%mask(d<=0) = 0; %mask(d>0) = 1; %image(uint8(mask .* I));colorbar;colormap(gray(256));

    mu_I = mean(I(:));
    std_I = std(I(:));

    %rows 740,680
    % 1220 1040
    %I = I(680:1160,740:1220);
    mpp_i_n = -1;
    MPP_ITERATIONS = -1;
    %% PARALLEL 
    %parfor Num_z = 1:300     

    %% Energy terms for AC
    alpha = 1;
    type = 0; 
    gaussian_kernel_dims = 3;

    g = ac_gradient_map(I, alpha, type,gaussian_kernel_dims*0, 0.5*eye(ndims(I)), 'same');

    mu_dark = mu_I - std_I;


    %% Energy Terms for MPP
    black_spots = zeros(size(I));
    black_spots(I < mu_I - 2*std_I) = 1;

    black_spots = double(bwareaopen(black_spots,30));
    black_spots(I > mu_I + std_I) = -1;

    g1 = ac_gradient_map(I, 0.005, type,gaussian_kernel_dims, 1*eye(ndims(I)), 'same');

    %% Contour Terms
    countour_weight = 0.5*90; 
    delta_t = 2; 
    n_iters = 350; 
    g_shape = {'circle','rectangle','circle_array','chess_box'};
    dims = size(I);

    %% MPP Terms

    list_of_contours = [];
    Overlap_Threshold = 0.1;


    list_of_contours = [];
    while(mpp_i_n <= MPP_ITERATIONS)

        %% Birth Step        
         r_gdc = rand(1,5);
        [phi,expansion_weight] = get_random_phi(r_gdc,mpp_i_n,dims);

        if(show_result)
            figure; imshow(I,[]); hold on;  
            c = contours(phi,[0,0]); 
            if isempty(c), break; end
            if exist('h','var') && all(ishandle(h)), delete(h); end; 
            h = zy_plot_contours(c,'linewidth',2); 

        end

        phi_s = ac_polymer_model((I-mu_dark)*14/255,phi, expansion_weight,countour_weight, g, ...
                                 delta_t, n_iters, show_result); 


        c = contours(phi_s,[0,0]); 
        [birth_step_contours] = get_list_of_contours(c);


        for cn=1:length(birth_step_contours)
            rx = birth_step_contours(cn).rx;
            ry = birth_step_contours(cn).ry;
            if(length(rx) < 20)
                birth_step_contours(cn).E = Inf;
                continue;
            end
            [E] = Snake_Energy(rx',ry',double(black_spots),g1);
            birth_step_contours(cn).E = E;
            birth_step_contours(cn).Vp = 0;
            if(E < Inf)
                birth_step_contours(cn).max_x = max(rx);
                birth_step_contours(cn).max_y = max(ry);
                birth_step_contours(cn).min_x = min(rx);
                birth_step_contours(cn).min_y = min(ry);
            end
        end

        %% Death Step
        %Sort Ellipses by likelyhood. Descending Order (less neg first).
        good_index = find([birth_step_contours.E] ~= Inf);
        birth_step_contours = birth_step_contours(good_index);


        % Merge with previous results
        if(~length(birth_step_contours))
            mpp_i_n = mpp_i_n + 1;
            continue;
        end

        list_of_contours = [list_of_contours,birth_step_contours];
        [~,sorted_by_E] =sort([list_of_contours.E]);
        list_of_contours = list_of_contours(sorted_by_E);

        mpp_i_n = mpp_i_n + 1;
        continue;
        number_of_contours = length(list_of_contours);


        if(number_of_contours == 0)
            mpp_i_n = mpp_i_n + 1;
            continue;
        end


        sn = 1;
        while(sn <= number_of_contours)    
           indexes_R_candidates = find_overlap_snakes(sn,list_of_contours);

           if(~isempty(indexes_R_candidates))
               overlap_R = find_overlaping_percent(I,list_of_contours(sn), list_of_contours(indexes_R_candidates), Overlap_Threshold);
           else
               overlap_R=0;
           end

           list_of_contours(sn).Vp = overlap_R;
           l = list_of_contours(sn).E;

           if(overlap_R > Overlap_Threshold)
               d_rate = 1.01;
           else
                d_rate = 0;
           end

           pp = rand(1,1);
           if(pp < d_rate)
               list_of_contours = list_of_contours([1:sn-1,sn+1:end]);
               number_of_contours = number_of_contours - 1;
           else
               sn = sn+1;
           end


        end

        mpp_i_n = mpp_i_n + 1;
    end


    %% Final Results

    im_segged = zeros(size(I));
    for i=1:length(list_of_contours)
        rx = list_of_contours(i).rx;
        ry = list_of_contours(i).ry;
        im_segged = DrawSegmentedArea2D(im_segged,rx,ry);
    end

    im_segged= im_segged';
    if(show_result)
        figure;
        [color_img] = plot_all_contours( I, list_of_contours,1);
        title('Final List');
    %else
        %[color_img] = plot_all_contours( I, list_of_contours);
    end

    if(Num_z > 999)
         output_name = ['SEGMENTED_VOIDS/voids_' num2str(Num_z) '.tiff'];
    elseif(Num_z > 99)
        output_name = ['SEGMENTED_VOIDS/voids_0' num2str(Num_z) '.tiff'];
    elseif(Num_z > 9)
        output_name = ['SEGMENTED_VOIDS/voids_00' num2str(Num_z) '.tiff'];
    else
        output_name = ['SEGMENTED_VOIDS/voids_000' num2str(Num_z) '.tiff'];
    end
    
    imwrite(im_segged, output_name);
    disp('Done');
    toc
end
 

