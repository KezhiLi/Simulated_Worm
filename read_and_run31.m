% save relative angles, central positions, and absolute angle

hdf5_file_name = 'Z:\Results\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11_skeletons.hdf5';
masked_image_file =  'Z:\MaskedVideos\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11.hdf5';

% read skeleton
ske_info = h5info(hdf5_file_name, '/skeleton');
frame_size = ske_info.Dataspace.Size(1:2);
frame_total = ske_info.Dataspace.Size(3);
ske_full = h5read(hdf5_file_name,'/skeleton'); 

% read stage position
stage_info = h5info(hdf5_file_name, '/stage_movement');
stage_vec = h5read(hdf5_file_name, '/stage_movement/stage_vec');
is_stage_move = h5read(hdf5_file_name, '/stage_movement/is_stage_move');
    pixelPerMicronX = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
    pixelPerMicronY = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_y');
    
stage_pixel = stage_vec(:,1:3:end);
stage_pixel(1,:) = abs(stage_pixel(1,:)/pixelPerMicronX);
stage_pixel(2,:) = abs(stage_pixel(2,:)/pixelPerMicronY);

ske = ske_full(:,1:3:end,1:3:end);
ske_diff = ske(:,2:end,:)-ske(:,1:end-1,:);


cut_ske = ske;
cut_stage = stage_pixel;

% fill the NaN with 1st order interpolation
cut_stage_fill = cut_stage;
mid_ske_ind = round((size(cut_ske,2)+1)/2);
for ii = 2: size(cut_stage_fill,2);
    if isnan(sum(cut_stage_fill(:,ii)))
        if ii == 2
            [non_val, non_ind] = min(isnan(cut_stage_fill(1,:)));
            cut_stage_fill(:,ii) = cut_stage_fill(:,non_ind);
        else
            %stage_mov = cut_ske(:,mid_ske_ind,ii)-2*cut_ske(:,mid_ske_ind,ii-1)+cut_ske(:,mid_ske_ind,ii-2);
            stage_mov = mean(cut_ske(:,:,ii),2)-2*mean(cut_ske(:,:,ii-1),2)+mean(cut_ske(:,:,ii-2),2);
            cut_stage_fill(:,ii) = (cut_stage_fill(:,ii-1) - stage_mov) + (cut_stage_fill(:,ii-1)-cut_stage_fill(:,ii-2));
        end
    end
end

% add global position to local position
sample_ske =zeros(size(cut_ske));
for ii = 1:size(cut_ske,3)
    sample_ske(:,:,ii) = cut_ske(:,:,ii) + kron(ones(1,size(cut_ske,2)),cut_stage_fill(:,ii));
end

% % repeat the vidoe 
% sample_ske_ori = sample_ske;
% for dd= 1:63;
%     sample_ske_next = sample_ske_ori;
%     % pixel compensation
%     %pixel_comp = kron(ones(1,size(sample_ske,2)),( -sample_ske(:,mid_ske_ind,1)+2*sample_ske(:,mid_ske_ind,end)-sample_ske(:,mid_ske_ind,end-1)   ));
%     pixel_comp = kron(ones(1,size(sample_ske,2)),( -mean(sample_ske(:,:,1),2)+2*mean(sample_ske(:,:,end),2)-mean(sample_ske(:,:,end-1),2) ));
%     for ii = 1:size(sample_ske_next,3);
%         sample_ske_next(:,:,ii) = sample_ske_next(:,:,ii) + pixel_comp;
%     end
%     sample_ske = cat(3, sample_ske, sample_ske_next);
% end
sample_ske_mid = squeeze(sample_ske(:,mid_ske_ind,:));
sample_ske_mid_diff = [[0;0],sample_ske_mid(:,2:end)-sample_ske_mid(:,1:end-1)];

% % show the animation of skeleton 
% for ii = 1:size(sample_ske,3);
%     %img = zeros(480,640);
%     plot(sample_ske(1,:,ii),sample_ske(2,:,ii),'*-');
%     axis equal
%     axis([2500 3600 6300 6600])
%     pause(0.2)
% end

sample_skediff =  ske_diff;

% convension from cart to polar coordinates
angle_vec = zeros(size(sample_skediff,2),size(sample_skediff,3));
len_vec = angle_vec;
for ii = 1:size(sample_skediff,3);
    [angle_vec(:,ii), len_vec(:,ii)] = cart2pol(sample_skediff(1,:,ii),sample_skediff(2,:,ii));
end

% calculate 1st row angle vector, or the 1st derivitive of the 1st row
% angle vector along time
firstRow_angle_vec = angle_vec(1,:);
firstRow_angle_vec_dif = [firstRow_angle_vec(1), firstRow_angle_vec(2:end) - firstRow_angle_vec(1:end-1)];
big_ind0 = firstRow_angle_vec_dif>pi;
firstRow_angle_vec_dif(big_ind0) = firstRow_angle_vec_dif(big_ind0)-2*pi;
sma_ind0 = firstRow_angle_vec_dif<-pi;
firstRow_angle_vec_dif(sma_ind0) = firstRow_angle_vec_dif(sma_ind0)+2*pi;    

% calculate the difference of angle vector along worm body
diff_angle_vec = angle_vec(2:end,:) - angle_vec(1:end-1,:);
big_ind = diff_angle_vec>pi;
diff_angle_vec(big_ind) = diff_angle_vec(big_ind)-2*pi;
sma_ind = diff_angle_vec<-pi;
diff_angle_vec(sma_ind) = diff_angle_vec(sma_ind)+2*pi;      


angle_vec3 = [diff_angle_vec;sample_ske_mid_diff;firstRow_angle_vec_dif];
%angle_vec3 = [diff_angle_vec;sample_ske_mid_diff;firstRow_angle_vec];

angle_vec3_ind = setdiff( 1:size(angle_vec3,2), find(isnan(sum(angle_vec3))==1));
angle_vec3 = angle_vec3(:,angle_vec3_ind);

%% generate data
% save('angle_vec.mat','angle_vec');
% save('angle_vec.xlsx','angle_vec');

% save_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
% if exist([save_path,'angle_vec3-2.hdf5'])>1
%     delete([save_path,'angle_vec3-2.hdf5']);
% end
% h5create([save_path,'angle_vec3-2.hdf5'], '/angle_vec3', size(angle_vec3), 'Datatype', 'double', ...
%     'Chunksize', size(angle_vec3), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
% h5write([save_path,'angle_vec3-2.hdf5'], '/angle_vec3', angle_vec3);

