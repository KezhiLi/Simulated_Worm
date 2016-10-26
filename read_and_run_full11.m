hdf5_file_name = 'Z:\Results\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11_skeletons.hdf5';
masked_image_file =  'Z:\MaskedVideos\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11.hdf5';

ske_info = h5info(hdf5_file_name, '/skeleton');
frame_size = ske_info.Dataspace.Size(1:2);
frame_total = ske_info.Dataspace.Size(3);
ske_full = h5read(hdf5_file_name,'/skeleton'); 

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
% 
% % for ii = 1:size(ske,3);
% %     img = zeros(480,640);
% %     plot(ske(1,:,ii),480-ske(2,:,ii),'*-');
% %     axis equal
% %     pause(0.5)
% % end
% 
% %% try to find the repeated motion that can use 
% % fra_try = 200;
% % cor_mtx = ones(fra_try,fra_try)*1e4;
% % for ii = 1:fra_try;
% %     for jj = ii+3:fra_try;
% %         cor_mtx(ii,jj)=sum(sum(abs(ske_diff(1,:,ii)-ske_diff(1,:,jj)).^2))...
% %             +sum(sum(abs(ske_diff(2,:,ii)-ske_diff(2,:,jj)).^2));
% %     end
% % end
% % 
% % figure, imshow(cor_mtx,[]);
% % 
% % cor_mtx_sm = cor_mtx.*(cor_mtx<30);
% % figure, imshow(cor_mtx_sm,[])
% % % similar frame index: 96, 140
% 
% cut_ske = ske(:,:,96:139);
% cut_stage = stage_pixel(:,96:139);
% sample_ske =zeros(size(cut_ske));
% for ii = 1:size(cut_ske,2)
%     sample_ske(:,:,ii) = cut_ske(:,:,ii) + kron(ones(1,size(cut_ske,2)),cut_stage(:,ii));
% end
% 
% sample_ske = cut_ske;
% for dd = 1:6;
%     sample_ske = cat(3, sample_ske, sample_ske);
% end
% % % show the animation of skeleton 
% % for ii = 1:size(sample_ske,3);
% %     img = zeros(480,640);
% %     plot(sample_ske(1,:,ii),480-sample_ske(2,:,ii),'*-');
% %     axis equal
% %     pause(0.2)
% % end

cut_skediff = ske_diff; %(:,:,96:139);
cut_skediff_ind = 1;
sample_skediff =[];
for ii = 1:size(cut_skediff ,3)
    if ~isnan(sum(sum(cut_skediff(:,:,ii))))
        sample_skediff(:,:,cut_skediff_ind) = cut_skediff(:,:,ii);
        cut_skediff_ind = cut_skediff_ind+1;
    end
end
%sample_skediff = cut_skediff;


% % repeat the process with 2^dd fold
% for dd = 1:7;
%     sample_skediff = cat(3, sample_skediff, sample_skediff);
% end

% convension from cart 2 polar coordinates
angle_vec = zeros(size(sample_skediff,2),size(sample_skediff,3));
len_vec = angle_vec;
for ii = 1:size(sample_skediff,3);
    [angle_vec(:,ii), len_vec(:,ii)] = cart2pol(sample_skediff(1,:,ii),sample_skediff(2,:,ii));
end
angle_vec_ori = angle_vec;

% make angle_vec continious
for ii = 1:size(angle_vec,1);
    for jj = 1:size(angle_vec,2)-1;
        while angle_vec(ii,jj)- angle_vec(ii,jj+1)>4 
            angle_vec(ii,jj+1) = 2*pi + angle_vec(ii,jj+1);
        end
        while  angle_vec(ii,jj)- angle_vec(ii,jj+1) < -4
            angle_vec(ii,jj+1) = -2*pi + angle_vec(ii,jj+1);
        end
    end
end

% transmitt the difference of angle_vec to ML
n_prev = 15;
angle_vec_diff = angle_vec(:,2:end) - angle_vec(:,1:end-1);
angle_vec_1c = angle_vec_ori(:,round(0.9*size(angle_vec,2)) + n_prev);

% %% debug use: show figure using components
% radias_vec=(A{n,1}*icasig{n,1});
% rho = median(len_vec,2);
%     pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
%     pred_ske = zeros(2,size(radias_vec,1)+1,size(radias_vec,2));
%     for ii = 1:size(radias_vec,2);
%         [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho);
%         pred_ske(:,2:end,ii) = cumsum(pred_ske_diff(:,:,ii),2);
%     end
%     for ii = 1:size(pred_ske,3);
%         %img = zeros(480,640);
%         plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
%         axis equal
%         pause(0.2)
%        % mov(ii) = save_crt_fra(filename_gif,ii, fps);
%     end
%     

% %% generate data
% % save('angle_vec.mat','angle_vec');
% % save('angle_vec.xlsx','angle_vec');
% % 
save_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
h5create([save_path,'angle_vec_full11.hdf5'], '/angle_vec_diff', size(angle_vec_diff), 'Datatype', 'double', ...
    'Chunksize', size(angle_vec_diff), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
h5write([save_path,'angle_vec_full11.hdf5'], '/angle_vec_diff', angle_vec_diff);
h5create([save_path,'angle_vec_full11.hdf5'], '/len_vec', size(len_vec), 'Datatype', 'double', ...
    'Chunksize', size(len_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
h5write([save_path,'angle_vec_full11.hdf5'], '/len_vec', len_vec);

h5create([save_path,'angle_vec_full11.hdf5'], '/angle_vec_1c', size(angle_vec_1c), 'Datatype', 'double', ...
    'Chunksize', size(angle_vec_1c), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
h5write([save_path,'angle_vec_full11.hdf5'], '/angle_vec_1c', angle_vec_1c);