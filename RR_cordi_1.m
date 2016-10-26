hdf5_file_name = 'Z:\Results\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11_skeletons.hdf5';
masked_image_file =  'Z:\MaskedVideos\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11.hdf5';

ske_info = h5info(hdf5_file_name, '/skeleton');
frame_size = ske_info.Dataspace.Size(1:2);
frame_total = ske_info.Dataspace.Size(3);
ske_full = double(h5read(hdf5_file_name,'/skeleton')); 

stage_info = h5info(hdf5_file_name, '/stage_movement');

stage_vec = h5read(hdf5_file_name, '/stage_movement/stage_vec');
is_stage_move = h5read(hdf5_file_name, '/stage_movement/is_stage_move');
    pixelPerMicronX = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
    pixelPerMicronY = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_y');
    
stage_pixel = stage_vec(:,1:3:end);
stage_pixel(1,:) = abs(stage_pixel(1,:)/pixelPerMicronX);
stage_pixel(2,:) = abs(stage_pixel(2,:)/pixelPerMicronY);

ske = ske_full(:,1:3:end,1:3:end);
ske_len = size(ske,2);
for ii = 1:size(stage_pixel,2)
    ske_adj(:,:,ii) =  ske(:,:,ii) +  stage_pixel(:,ii)*ones(1,ske_len);
end

ske_NoNan = ske_adj(:,:,~isnan(ske_adj(1,1,:)));
min_x = min(min(ske_NoNan(1,:,:)));
min_y = min(min(ske_NoNan(2,:,:)));
ske_NoNan_adj(1,:,:) = ske_NoNan_adj(1,:,:)- min_x;
ske_NoNan_adj(2,:,:) = ske_NoNan_adj(2,:,:)- min_y;
xy_adj = [min_x ,min_y];

ske_vec =reshape(ske_NoNan, size(ske_NoNan,1)*size(ske_NoNan,2),size(ske_NoNan,3));

cord_1 = [ske_vec;ones(1,size(ske_vec,2))];

save_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
if exist([save_path,'cord_1.hdf5'])>1
    delete([save_path,'cord_1.hdf5']);
end

h5create([save_path,'cord_1.hdf5'], '/cord_1', size(cord_1), 'Datatype', 'double', ...
    'Chunksize', size(cord_1), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
h5write([save_path,'cord_1.hdf5'], '/cord_1', cord_1);

h5create([save_path,'cord_1.hdf5'], '/xy_adj', size(xy_adj), 'Datatype', 'double', ...
    'Chunksize', size(xy_adj), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
h5write([save_path,'cord_1.hdf5'], '/xy_adj', xy_adj);

% % show the animation of skeleton 
% for ii = 1:size(ske_NoNan,3);
%     %img = zeros(480,640);
%     plot(ske_NoNan(1,:,ii),ske_NoNan(2,:,ii),'*-');
%     axis equal
%     pause(0.2)
% end