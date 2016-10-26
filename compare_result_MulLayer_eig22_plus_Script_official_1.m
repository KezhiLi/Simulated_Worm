% 6 dimensions + 1 DiffOfMeanofAbsAngles


addpath('X:\Kezhi\fastICA');


fps = 5;
% 
% csv_predicted = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_predicted_full22.csv';
% %csv_predicted = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\247 JU438 on food L_2011_03_03__11_18___3___1_predicted-hid4.csv';
% predicted = csvread(csv_predicted);
% csv_y_test = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_y_test_full22.csv';
% %csv_y_test = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\04-03-11\2_y_test_rev2.csv';
% y_test = csvread(csv_y_test);
% %csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_generated_full22.csv';
% %csv_generated = 'Z:\DLWeights_temp\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\247 JU438 on food L_2011_03_03__11_18___3___1_generated-hid4.csv';
% %csv_generated = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-18\!worm_videos\copied from pc207-13\misc_videos\robyn dvd\AQ2000 on food L_2011_05_19__16_42_00___1___12_eig(4hid).csv' ;

csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\test_res_2209_2016\muiltiFile_600epoch\generate_round5\1_generated_long1.csv';
%csv_generated = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\1_relu_generated_1.csv';
generated_ske = csvread(csv_generated);

% load_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
% eig_vec_file = 'eig_para_full22.hdf5';
load_path = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\';
eig_vec_file = '247 JU438 on food L_2011_03_03__11_18___3___1_eig.hdf5'; %1 
%eig_vec_file = '431 JU298 on food L_2011_03_03__12_29___3___5_eig.hdf5';%2
%eig_vec_file = '764 ED3049 on food L_2011_03_03__15_44___3___7_eig.hdf5';%3 
%eig_vec_file = '800 MY16 on food L_2011_03_03__16_52___3___11_eig.hdf5';%4
%eig_vec_file = '972 JU345 on food L_2011_03_03__11_53___3___3_eig.hdf5';%5
%eig_vec_file = 'N2 on food R_2011_03_03__16_17___3___9_eig.hdf5';%6

eig_vec = h5read([load_path,eig_vec_file],'/eig_vec'); 
%eig_vec = A{n};
len_vec = h5read([load_path,eig_vec_file],'/len_vec'); 
mean_angle_vec = h5read([load_path,eig_vec_file],'/mean_angle_vec'); 

FirstNoFrm = round(size(len_vec,2)*0.9+50);
FirstAbsAng = mean_angle_vec(FirstNoFrm);

for tt = 3;
    if tt==1
        eig_radias_vec = predicted(:,1:end-1);
        mean_angle_vec_diff = predicted(:,end);
        mean_angle_vec_cur = (mean_angle_vec(FirstNoFrm:end-1))'+mean_angle_vec_diff;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_1(22_4hid).gif';
         
        eig_radias_vec2 = y_test(:,1:end-1);
        mean_angle_vec_diff2 = y_test(:,end);
        mean_angle_vec_diff2(1) = mean_angle_vec_diff2(1) + FirstAbsAng;
        mean_angle_vec_cur2 = cumsum(mean_angle_vec_diff2);

    elseif tt ==2
    %    eig_radias_vec = icasig{n};
        eig_radias_vec = y_test(:,1:end-1);
        mean_angle_vec_diff = y_test(:,end);
        mean_angle_vec_diff(1) = mean_angle_vec_diff(1) + FirstAbsAng;
        mean_angle_vec_cur = cumsum(mean_angle_vec_diff);
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\muil_eig_ytest(rev2-0)_F2.gif';
        
        eig_radias_vec = eig_radias_vec(1:400,:);
        mean_angle_vec_cur =mean_angle_vec_cur(1:400,:);
    elseif tt == 3
        eig_radias_vec = generated_ske(:,1:end-1);
        mean_angle_vec_diff = generated_ske(:,end);
        mean_angle_vec_diff(1) = mean_angle_vec_diff(1) + FirstAbsAng;
        mean_angle_vec_cur = cumsum(mean_angle_vec_diff);
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\200_200_200_200_600ep_6.gif';
    end
   
    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end
    
    radias_vec=eig_vec*eig_radias_vec+ (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur))';

    % radias to ske
    rho = median(len_vec,2);
    pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
    pred_ske = zeros(2,size(radias_vec,1)+1,size(radias_vec,2));
    for ii = 1:size(radias_vec,2);
        [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho);
        pred_ske(:,2:end,ii) = cumsum(pred_ske_diff(:,:,ii),2);
        oriPoint = floor(size(pred_ske,2)+1)/2;
      %  pred_ske(:,:,ii) = pred_ske(:,:,ii)- kron(ones(1,size(pred_ske,2)),pred_ske(:,oriPoint,ii));
      pred_ske(:,:,ii) = pred_ske(:,:,ii)- kron(ones(1,size(pred_ske,2)),mean(pred_ske(:,:,ii),2));
    end
end


% please change the first load() to some file saved in local. 
% you can switch between 'load('pred_ske_1.mat');' or
% 'load('pred_ske_1.mat')', to load skeleton x,y coodidates without or with
% mean angles in consideration

addpath('X:\Andre\eric\RabetsEtAlModel - Copy\');

%% read data
% % load arbitrary worm data
load(['Z:/single_worm/samples/switched_sample/' ...
    'acr-21 (ok1314)III on food L_2010_02_24__14_45_13__11_features.mat'])
  %  'osm-9 (ky10) on food R_2010_06_15__14_57_24___8___8_features.mat'])

% % read pred_ske without mean angles 
% load('pred_ske_4.mat');

% % read pred_ske with mean angles
% load('pred_ske_3.mat'); % or 'pred_ske_3.mat' or 'pred_ske_4.mat'
% ------------------------------Part A-------------------------------------

% plot two skeletons to show motion
frame1 = 20;  % 1205
frame2 = 25;  % 1230

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

% 
% figure
% 
% div = 4;
% Xrecon_short = Xrecon(1:round(size(Xrecon,1)/div),:);
% Yrecon_short = Yrecon(1:round(size(Xrecon,1)/div),:);
% %% show the animation of skeleton
% for ii = 1:size(Xrecon_short,1);
%     
%     plot(Xrecon_short(ii, :)' - Xrecon_short(1, 1), Yrecon_short(ii, :)' - Yrecon_short(1, 1), 'Color', 'r')
%     hold on
%     plot(Xrecon_short(ii, 1)' - Xrecon_short(1, 1), Yrecon_short(ii, 1)' - Yrecon_short(1, 1), ...
%         '.', 'Color', 'b', 'MarkerSize', 15)
%     axis equal
%     xlim([ min(min(Xrecon_short - Xrecon_short(1, 1))), max(max(Xrecon_short - Xrecon_short(1, 1)))])
%     ylim([ min(min(Yrecon_short - Yrecon_short(1, 1))), max(max(Yrecon_short - Yrecon_short(1, 1)))])
%     
%     hold off
%     pause(0.05)
%       %    mov(ii) = save_crt_fra(filename_gif,ii, fps);
% end

figure, 
plot(Xrecon(:, :)' - Xrecon(1, 1), Yrecon(:, :)' - Yrecon(1, 1), 'Color', 'r')
hold on
%plot(Xrecon(:, 1)' - Xrecon(1, 1), Yrecon(:, 1)' - Yrecon(1, 1), ...
%    '.', 'Color', 'b', 'MarkerSize', 15)
axis equal
xlim([ min(min(Xrecon - Xrecon(1, 1))), max(max(Xrecon - Xrecon(1, 1)))])
ylim([ min(min(Yrecon - Yrecon(1, 1))), max(max(Yrecon - Yrecon(1, 1)))])

hold off


% saveas(gcf,'C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\08-03-11\1000epochs\large_noise\N2 on food L_2011_03_08__15_53___3___9_eig.fig')