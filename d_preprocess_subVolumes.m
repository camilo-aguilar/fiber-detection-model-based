%% File: preprocess_subVolumes
%  Author: Camilo Aguilar
%  Funciton: Assign a number to each detected fiber and store it
%
%  Input:  Folders ../SUBVOLUMES/sVn (where n = 1,2,3,...75)
%          Each sVn folder must have:
%               fibers/mer_i.tiff (450 images) (i=1,...,450)
%               fibers_info/merged_fiber.txt
%
%  Output: fibers_info/Vn.mat with:
%               1)Vn :      450x450x450 double with pixels denoting fiber # 
%               2)Vn_info:  text file containing the fiber # and info
%               
%% Pre-Process Subvolumes and their info
addpath('include');
tic
offset = 0;
num_fibers = 0;
for SubVn=1:75 
    disp(SubVn);
    % Read volumes and text file:
    txt_dir = ['SUBVOLUMES/sV' num2str(SubVn) '/fibers_info/merged_fiber.txt'];
    img_dir = ['SUBVOLUMES/sV' num2str(SubVn) '/fibers/mer_'];
    V = get_subVolume( img_dir, 450);
    fileID = fopen(txt_dir,'r');
    Vn_info = fscanf(fileID,'%d,%f,%f,%f,%f,%f,%f,%f',[8 Inf]);
    Vn_info = Vn_info(:,2:end);
    Vn_info(1,:) = Vn_info(1,:) - 1;
    fclose(fileID);
    
    % Transform volume from rgb to double
    Vn = squeeze(V(:,:,3,:) + 256 * V(:,:,2,:) + (256^2 * V(:,:,1,:)));
    num_fibers = Vn_info(1,end);
    num_fibers2 = max(Vn(:));
    if(num_fibers2 ~= num_fibers)
        disp(['Error in Volume ' num2str(SubVn)]);
        return;
    end
    
    % Give Offset from previous volumes
    Vn = Vn + offset .* sign(Vn);
    Vn_info(1,:) = Vn_info(1,:) + offset;
    offset = offset + num_fibers;
    processed_dir = ['SUBVOLUMES/sV' num2str(SubVn) '/fibers_info/'];
    save([processed_dir 'Vn.mat'],'Vn','Vn_info');
    
end
toc