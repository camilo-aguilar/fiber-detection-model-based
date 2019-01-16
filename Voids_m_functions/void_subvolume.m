function[ ] = void_subvolume( SubVn )
%% function: void_subvolume(SubVn)
    %  Author: Camilo Aguilar
    %  Funciton: Find voids in the subvolume number: SubVn
    %
    %  Parameter: Subvn: volume number
    %
    %  Required:    Folders ../SUBVOLUMES/sVn (where n = 1,2,3,...75)
    %               Each sVn folder must have:
    %                       /data/subV_i.tiff (450 images) (i=1,...,450)
    %                       /seg/subV_seg_i.tiff (450 images) (i=1,...,450)
    %
    %  Output: 450 images with: ../SUBVOLUMES/sVn/
    %                           voids/subVn_void_i (i=1,...,450)
    %% Pre-Process Subvolumes and their info
    addpath('util_functions');
    
    num_images = 450;
    
    show_result = 0;
    alpha = 1;
    type = 0; 
    gaussian_kernel_dims = 3;
    MPP_ITERATIONS = 2;

    %% Contour Terms
    countour_weight = 0.5*90; 
    delta_t = 2; 
    n_iters = 600; 
    min_length = 20;
    
    %% MPP Terms
    % Overlap_Threshold = 0.1;
    
    %% Input Images    
    dir_v = ['../SUBVOLUMES/sV' num2str(SubVn)];
    if(~exist([dir_v '/voids'], 'dir'))
        mkdir([dir_v '/voids']);
    end
    
    volume_voids = zeros(450,450,num_images);
    for Num_z=1:num_images
        disp(Num_z);
        % Read Images
        if(Num_z) <10
            I = imread([dir_v '/data/subV_000' num2str(Num_z) '.tif']);
            Im_seg = imread([dir_v '/seg/subV_seg_000' num2str(Num_z) '.tif']);
        elseif(Num_z) <100
            I = imread([dir_v '/data/subV_00' num2str(Num_z) '.tif']);
            Im_seg = imread([dir_v '/seg/subV_seg_00' num2str(Num_z) '.tif']);
        elseif(Num_z) <1000
            I = imread([dir_v '/data/subV_0' num2str(Num_z) '.tif']);
            Im_seg = imread([dir_v '/seg/subV_seg_0' num2str(Num_z) '.tif']);
        else
            I = imread([dir_v '/data/subV_' num2str(Num_z) '.tif']);
            Im_seg = imread([dir_v '/seg/subV_seg_' num2str(Num_z) '.tif']);
        end
        
        I = double(I);
        [rows,cols] = size(I);

        mu_I = mean(I(:));
        std_I = std(I(:));


        %% Energy terms for AC
        g = ac_gradient_map(I, alpha, type,gaussian_kernel_dims*0, 0.5*eye(ndims(I)), 'same');
        mu_dark = mean(I(Im_seg < 90));
        
        %% Energy Terms for MPP
        black_spots = zeros(size(I));
        black_spots(Im_seg == round(255/6)) = 1;

        black_spots = double(bwareaopen(black_spots,40));
        black_spots(Im_seg > 3* round(255/6)) = -1;

        g1 = ac_gradient_map(I, 0.005, type,gaussian_kernel_dims, 1*eye(ndims(I)), 'same');

        dims = size(I);  

        list_of_contours = [];
        mpp_i_n = 1;
        while(mpp_i_n <= MPP_ITERATIONS)
            %% Birth Step        
             r_gdc = rand(1,5);
            [phi,expansion_weight] = get_random_phi(r_gdc,mpp_i_n,dims);

            phi_s = ac_polymer_model((I-mu_dark)*14/255,phi, expansion_weight,countour_weight, g, ...
                                     delta_t, n_iters, show_result); 

            c = contours(phi_s,[0,0]); 
            [birth_step_contours] = get_list_of_contours(c);

            for cn=1:length(birth_step_contours)
                rx = birth_step_contours(cn).rx;
                ry = birth_step_contours(cn).ry;
                
                if(length(rx) < min_length)
                    birth_step_contours(cn).E = Inf;
                    continue;
                end
                [E] = Snake_Energy(rx',ry',double(black_spots),g1);
                birth_step_contours(cn).E = E;
                
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

            
            
            %%%%%%% Death Step is skipped %%%%%%%%
            mpp_i_n = mpp_i_n + 1;
            continue;
            %%%%%%% Delete this part if keep Death Step  %%%%%%%%

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
        im_segged = min(1, im_segged);

        volume_voids(:,:,Num_z) = im_segged;
    end

    V = imgaussfilt3(volume_voids - 0.5,3);

    V = uint8(255 .* sign(V));

    for i=1:num_images
         if(i) <10
            name = [dir_v '/voids/subVn_void_000' num2str(i) '.tiff'];
        elseif(i) <100
            name = [dir_v '/voids/subVn_void_00' num2str(i) '.tiff'];
        elseif(i) <1000
            name = [dir_v '/voids/subVn_void_0' num2str(i) '.tiff'];
        else
            name = [dir_v '/voids/subVn_void_' num2str(i) '.tiff'];
        end
        imwrite(uint8(V(:,:,i)), name);
    end
    

end