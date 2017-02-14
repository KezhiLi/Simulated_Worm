% run Compare_models_1.m to run python scripts given RNN model, generate
% pridction .csv, and draw the trajectories based on worm 'shape to
% trajectory' function, then save the .fig and .png trajectory figure to
% local folder
%
% USE newer version codes from Eric to convert posture to 2-d locomotion
%
% Kezhi Li, 13st Feb, 2016
%

%%please change stain here, and in the python file 'Function_predict_multiLayer_GiveModel_input.py'
stain = 'N2';    %

%% run the algorithm
%cur_folder = 'Z:\DLWeights\test_files_folder\';
model_path = ['C:/Users/kezhili/Documents/Python Scripts/data/FromAWS/',stain,'/multiFile_',stain,'_7-260-260-260-260-7_600ep.h5'];
%cur_folder_input = 'Z:/DLWeights/eig_catagory_Straits/ED3049/';
cur_folder = ['Z:\DLWeights\eig_catagory_Straits\',stain,'\'];
% 6 dimensions + 1 DiffOfMeanofAbsAngles

% %use python
% %commandStr = ['python "C:\Users\kezhili\Documents\Python Scripts\Function_predict_multiLayer_GiveModel_input.py" ', '"',model_path,'"', '"',cur_folder,'"'];
% commandStr = ['python "C:\Users\kezhili\Documents\Python Scripts\Function_predict_multiLayer_GiveModel_input_NewVer.py" ', '"',model_path,'"'];
% system(commandStr)
% disp('python done!')



addpath('X:\Kezhi\fastICA');
addpath('C:\Kezhi\MyCode!!!\Simulated_Worm\Test_Eric\');

fps = 5;

% find all .fig file names
all_csv_file = subdir([cur_folder,'*.csv']);

num_csv = size(all_csv_file,1);

for nf = 1:num_csv;  % 476
    disp([num2str(nf),'/',num2str(num_csv)])
    csv_file_name = all_csv_file(nf).name

    csv_generated = csv_file_name;
    generated_ske = csvread(csv_generated);

    % read an arbitrary file
    load_path = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\';
    eig_vec_file = '247 JU438 on food L_2011_03_03__11_18___3___1_eig.hdf5'; %1 

    eig_vec = h5read([load_path,eig_vec_file],'/eig_vec'); 
    len_vec = h5read([load_path,eig_vec_file],'/len_vec'); 
    mean_angle_vec = h5read([load_path,eig_vec_file],'/mean_angle_vec'); 

    FirstNoFrm = round(size(len_vec,2)*0.9+50);%length(mean_angle_vec)-size(y_test,1);
    FirstAbsAng = mean_angle_vec(FirstNoFrm);

    eig_radias_vec = generated_ske(:,1:end-1);
    mean_angle_vec_diff = generated_ske(:,end);
    mean_angle_vec_diff(1) = mean_angle_vec_diff(1) + FirstAbsAng;
    mean_angle_vec_cur = cumsum(mean_angle_vec_diff);

    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end

    radias_vec=eig_vec*eig_radias_vec; % + (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur))';

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
    [DX, DY, ODX, ODY, VX, VY, Xtil, Ytil, THETA] = ...
        subtractRBM(X, Y, XCM, YCM, UX, UY, UXCM, UYCM, OMEG, 1);
    
    % rotate velocites and tangent angles to the worm frame of reference
    [TX, TY] = lab2body(TX, TY, THETA);
    [VX, VY] = lab2body(VX, VY, THETA);
    
    % use model to predict rigid body motion
    alpha = 100; % large alpha corresponds to no-slip during crawling
    RBM = posture2RBM(TX, TY, Xtil, Ytil, VX, VY, 1, I, ds, alpha);
    
    % calculate the predicted rigid body motion
    [XCMrecon, YCMrecon, THETArecon] = integrateRBM(RBM, 1, THETA);
    
    % add the rigid body motion to the skeleton coordinates
    [Xrecon, Yrecon] = addRBMRotMat(Xtil, Ytil, XCMrecon, YCMrecon, ...
        THETArecon, XCM, YCM, THETA);
    
        % plot a movie of worm with motion back in
    for ii = 1:1:size(Xrecon,1) % only plotting every 10 frames to speed up movie
        plot(Xrecon(ii, :)- Xrecon(1, 1), Yrecon(ii, :) - Yrecon(1, 1), 'Color', 'r')
        hold on
        plot(Xrecon(ii,1)- Xrecon(1, 1), Yrecon(ii, 1) - Yrecon(1, 1), '.', 'Color', 'b', 'MarkerSize', 15)
        axis equal
        xlim([ min(min(Xrecon - Xrecon(1, 1))), max(max(Xrecon - Xrecon(1, 1)))])
        ylim([ min(min(Yrecon - Yrecon(1, 1))), max(max(Yrecon - Yrecon(1, 1)))])

        hold off
        pause(0.01) % forces refresh so you can see the plot
    end

%     figure, 
%     plot(Xrecon(:, :)' - Xrecon(1, 1), Yrecon(:, :)' - Yrecon(1, 1), 'Color', 'r')
%     hold on
% 
%     axis equal
%     xlim([ min(min(Xrecon - Xrecon(1, 1))), max(max(Xrecon - Xrecon(1, 1)))])
%     ylim([ min(min(Yrecon - Yrecon(1, 1))), max(max(Yrecon - Yrecon(1, 1)))])
% 
%     hold off
%     % save fig
%     saveas(gcf,[csv_file_name(1:end-3),'fig'])
%     % save png
%     print('-dpng', '-r600', strcat(csv_file_name(1:end-3),'png'));
%     close all

end

