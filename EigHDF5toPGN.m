% Convert _eig.hdf5 to .pgn of trajectories directly
%
% Kezhi Li, 18, Jan, 2017
%


%% run the algorithm

%cur_folder_input = 'Z:/DLWeights/eig_catagory_Straits/ED3049/';
cur_folder = 'C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\wild-isolate\';
% 6 dimensions + 1 DiffOfMeanofAbsAngles

addpath('X:\Kezhi\fastICA');
addpath('X:\Andre\eric\RabetsEtAlModel - Copy\');

fps = 5;

% find all .fig file names
all_hdf5_file = subdir([cur_folder,'*_eig.hdf5']);

num_hdf5 = size(all_hdf5_file,1);

for nf = 1:num_hdf5;  % 476
    disp([num2str(nf),'/',num2str(num_hdf5)])
    hdf5_file_name = all_hdf5_file(nf).name
    
    % read information
    ske_info = h5info(hdf5_file_name);
    
    eig_ske = h5read(hdf5_file_name,'/eig_coef'); 

    eig_vec = h5read(hdf5_file_name,'/eig_vec'); 
    len_vec = h5read(hdf5_file_name,'/len_vec'); 
    mean_angle_vec = h5read(hdf5_file_name,'/mean_angle_vec');

    eig_radias_vec = eig_ske(1:end-1,:);
    mean_angle_vec_diff = eig_ske(end,:);
    mean_angle_vec_cur = cumsum(mean_angle_vec_diff);

    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end

    radias_vec=eig_vec*eig_radias_vec+ (kron(ones(size(eig_vec,1),1), mean_angle_vec_cur));

    % radias to ske
    rho = median(len_vec,2);
    pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
    pred_ske = zeros(2,size(radias_vec,1)+1,size(radias_vec,2));
    for ii = 1:size(radias_vec,2);
        [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho);
        pred_ske(:,2:end,ii) = cumsum(pred_ske_diff(:,:,ii),2);
        oriPoint = floor(size(pred_ske,2)+1)/2;
        pred_ske(:,:,ii) = pred_ske(:,:,ii)- kron(ones(1,size(pred_ske,2)),mean(pred_ske(:,:,ii),2));
    end

    %% read data
    % % load arbitrary worm data
    load(['Z:/single_worm/samples/switched_sample/' ...
        'acr-21 (ok1314)III on food L_2010_02_24__14_45_13__11_features.mat'])
      %  'osm-9 (ky10) on food R_2010_06_15__14_57_24___8___8_features.mat'])

    % ------------------------------Part A-------------------------------------

    worm.posture.skeleton.x = squeeze(pred_ske(1,:,:));
    worm.posture.skeleton.y = squeeze(pred_ske(2,:,:));

    %worm.morphology.length(:) = mean(sum(len_vec));
    worm.morphology.length(:) = 1;

    % ------------------------------Part B-------------------------------------

    % plot series of skeletons, subtract away rigid body motion, plot model
    % result

    % first pre-process skeleton data
    [X, Y] = preprocSkel(worm.posture.skeleton.x, ...
        worm.posture.skeleton.y, nanmean(worm.morphology.length));

    % calculate arclength increment
    ds = 1/(size(X, 2)-1);

    % calculate the observed rigid body motion
    [XCM, YCM, UX, UY, UXCM, UYCM, TX, TY, NX, NY, I, OMEG] = ...
        getRBM(X, Y, 1, ds, 1);

    % subtract the observed rigid body motion
    [DX, DY, ODX, ODY, VX, VY] = ...
        subtractRBM(X, Y, XCM, YCM, UX, UY, UXCM, UYCM, OMEG);

    % use model to predict rigid body motion
    alpha = 86.1;
    RBM = posture2RBM(TX, TY, DX, DY, VX, VY, 1, I, ds, alpha);

    % add the predicted rigid body motion back to the shifted skeletons
    [VXmod, VYmod, Xrecon, Yrecon] = addRBMRot(DX, DY, VX, VY, RBM, 1);

    figure, 
    plot(Xrecon(:, :)' - Xrecon(1, 1), Yrecon(:, :)' - Yrecon(1, 1), 'Color', 'r')
    hold on

    axis equal
    xlim([ min(min(Xrecon - Xrecon(1, 1))), max(max(Xrecon - Xrecon(1, 1)))])
    ylim([ min(min(Yrecon - Yrecon(1, 1))), max(max(Yrecon - Yrecon(1, 1)))])

    hold off
%     % save fig
%     saveas(gcf,[hdf5_file_name(1:end-3),'fig'])
%     % save png
    print('-dpng', '-r600', strcat(hdf5_file_name(1:end-4),'png'));
    close all

end

