%function  ske2basicworm(cur_folder)

% 
% 
% 
% 
% 
% 
% 21st Oct 2016, Kezhi Li, Copyright@Imperial College

% cur_folder = 'Z:\DLWeights\test_files_folder\';
cur_folder = 'C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\7-200-200-200-200-7\07-03-11\600epochs\no_noise\';

addpath('X:\Kezhi\fastICA');
addpath('X:\Andre\eric\RabetsEtAlModel - Copy\');
addpath(genpath('C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_31\'));

% find all .csv file names
all_csv_file = subdir([cur_folder,'*.csv']);
num_csv = size(all_csv_file,1);

for nf = 1:num_csv;  
    disp([num2str(nf),'/',num2str(num_csv)])
    csv_file_name = all_csv_file(nf).name
    [pathstr,name1,ext] = fileparts(csv_file_name);

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

    radias_vec=eig_vec*eig_radias_vec+ (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur))';

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

    

    %% Calculate the skeleton
    no_frame = size(Xrecon,1);
    % the length of points on skeleton
    m_fre_pt = size(Xrecon,2);
    % vector prepared for skeleton
    t = 1:m_fre_pt ;
    seg_len = 8;
    ts = 1:1/(seg_len*2):m_fre_pt;
    width = 10;
    
    all_skeletons = cell(1,no_frame);
    all_non_vulva_contours = cell(1,no_frame);
    all_vulva_contours = cell(1,no_frame);
        
    for nn = 1: no_frame;

        worm_ske1 = [Xrecon(nn,:)',Yrecon(nn,:)'];
    %     % The skeleton curve, which is a combination of spline curve and linear curve
    %     worm_ske_curv_x =  0.7*interp1(t,worm_ske1(:,1),ts,'spline')+0.3*interp1(t,worm_ske1(:,1),ts,'linear') ;
    %     worm_ske_curv_y =  0.7*interp1(t,worm_ske1(:,2),ts,'spline')+0.3*interp1(t,worm_ske1(:,2),ts,'linear') ;

        ske_pred_xy1 = nakeinterp1(t',Xrecon(nn,:),ts');
        ske_pred_xy2 = nakeinterp1(t',Yrecon(nn,:),ts');

   %% Based on the skeleton curve and width, calculate the contour on both sides of the worm. 

        % vector prepared for contour   
        tc = 1:2*m_fre_pt ;
        tsc = 1:1/round(seg_len*2.3):2*m_fre_pt;
        
        % Calculate the perpendicular direction of each points on skeleton
        [worm_ske1_T,worm_ske1_N] = frenet_TN(Xrecon(nn,:)',Yrecon(nn,:)');
        % Calculate points on the contour
        [worm_shape,worm_body] = ske2shape(worm_ske1, worm_ske1_N, width, -0.4);
        % Calculate the curve links all the points on contour
        worm_shape_x =  nakeinterp1(tc',worm_body(:,1),tsc');
        worm_shape_y =  nakeinterp1(tc',worm_body(:,2),tsc');
        len_contour_tot = length(worm_shape_x);    
        %seperate contour to 'vulva_contours' and 'non_vulva_contours'
        if ~~(len_contour_tot/2-floor(len_contour_tot/2)) 
            half_ind = round(len_contour_tot/2);
        else
            half_ind = round(len_contour_tot/2)+1;
        end
%         %% draw
%         plot(ske_pred_xy1(1:end),ske_pred_xy2(1:end),'LineWidth',2,'Color',[1 0.1 0.1]);
%         hold on 
%         plot(worm_shape_x(1:end),worm_shape_y(1:end),'LineWidth',2,'Color',[0.1 1 0.1]);

        %% set result variables
        all_skeletons{nn} = single([ske_pred_xy1,ske_pred_xy2]);
        all_non_vulva_contours{nn} = single([worm_shape_x(1:half_ind),worm_shape_y(1:half_ind)]);
        all_vulva_contours{nn} = flip(single([worm_shape_x(half_ind:end),worm_shape_y(half_ind:end)]));
      
% calculate is_valid based on all_skeleton
        
    end
    is_valid = ~cellfun(@isempty,all_skeletons);
    is_stage_movement = [~is_valid,logical(0)];
    %% save info.mat
   save([pathstr,'\',name1,'_info.mat'], 'all_skeletons', 'all_non_vulva_contours', ...
       'all_vulva_contours', 'is_stage_movement', 'is_valid','-v7.3')
   

end


