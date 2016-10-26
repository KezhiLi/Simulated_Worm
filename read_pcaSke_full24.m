% read_pcaSke_full21
% 5 eigen coefficients
% 1 difference of mean of the absolute angle (after adjusted)

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

cut_skediff = ske_diff;
cut_skediff_ind = 1;
sample_skediff =[];
for ii = 1:size(cut_skediff ,3)
    if ~isnan(sum(sum(cut_skediff(:,:,ii))))
        sample_skediff(:,:,cut_skediff_ind) = cut_skediff(:,:,ii);
        cut_skediff_ind = cut_skediff_ind+1;
    end
end

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

% mean over all angles for a time
mean_angle_vec = mean(angle_vec);
diff_mean_angle_vec =[mean_angle_vec(1) ,mean_angle_vec(2:end) - mean_angle_vec(1:end-1)];


angle_vec_adj = angle_vec - kron(ones(size(angle_vec,1),1),mean_angle_vec);

%I am going to do different number of eigenvalues though.
n=5; %Number of eigenvalues. It doesn't do just that, it does them all up
     %till then.
icasig=cell(n,1); %These will be the actual values.
A=cell(n,1); %These are the eigenshapes.
W=cell(n,1); %These are the inverses of the eigenshapes matrix, they aren't 
             %actually interesting.
lasteig=1:n;



%DATAsrc = wildtypedata;
for i=1:n
    [icasig{i,1}, A{i,1}, W{i,1}]=fastica(angle_vec_adj,'approach','defl', ...
        'numOfIC', lasteig(1,i), 'g', 'pow3', 'finetune', 'pow3', ...
        'stabilization', 'on', 'lasteig', lasteig(1,i));
end

%give weights
icasig_n = icasig{n};
icasig_n(1,:) = icasig_n(1,:)*2;
icasig_n(2,:) = icasig_n(2,:)*2;


eig_coef = [icasig_n;diff_mean_angle_vec];
eig_vec = A{n};

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
% save_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
% h5create([save_path,'eig_para_full24.hdf5'], '/eig_coef', size(eig_coef), 'Datatype', 'double', ...
%     'Chunksize', size(eig_coef), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
% h5write([save_path,'eig_para_full24.hdf5'], '/eig_coef', eig_coef);
% h5create([save_path,'eig_para_full24.hdf5'], '/eig_vec', size(eig_vec), 'Datatype', 'double', ...
%     'Chunksize', size(eig_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
% h5write([save_path,'eig_para_full24.hdf5'], '/eig_vec', eig_vec);
% h5create([save_path,'eig_para_full24.hdf5'], '/len_vec', size(len_vec), 'Datatype', 'double', ...
%     'Chunksize', size(len_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
% h5write([save_path,'eig_para_full24.hdf5'], '/len_vec', len_vec);

