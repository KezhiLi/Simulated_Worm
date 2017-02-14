% function: show intermediam c matrices for a worm video 
% 
% 
% 
% 
% 

% strain name
stain =  'N2';
cur_folder = ['Z:\DLWeights\eig_catagory_Straits\',stain,'\'];
% the hdf5 file storages the c matrices
name = 'interm_N2 on food L_2011_02_17__12_51_07___7___7_eig_7-260-260-260-260-7_600ep.hdf5';
feature_name = 'N2 on food L_2011_02_17__12_51_07___7___7_features.hdf5';
hdf5_file_name = [cur_folder,name];


idx_layer = 1;
idx_neuron = 1;
% the index of layer and neuron one wants to observe
idx = [idx_layer, idx_neuron];

% c matrix of different layers
c_L = cell(4);
y_test = h5read(hdf5_file_name,'/y_test'); 
c_L{1} = h5read(hdf5_file_name,'/c_L0_col'); 
c_L{2} = h5read(hdf5_file_name,'/c_L1_col'); 
c_L{3} = h5read(hdf5_file_name,'/c_L3_col'); 
c_L{4} = h5read(hdf5_file_name,'/c_L5_col'); 
eig_vec = h5read(hdf5_file_name,'/eig_vec'); 
len_vec = h5read(hdf5_file_name,'/len_vec'); 

%c_L = [c_L0_col;c_L1_col;c_L3_col;c_L5_col];


addpath('X:\Kezhi\fastICA');
addpath('C:\Kezhi\MyCode!!!\Simulated_Worm\Test_Eric\');

for nf = 1  %num_csv;  % 476
    disp(name)
    % read eigen vectors (real worms)
    eig_radias_vec = y_test(1:end-1,:);

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
    
    %% plot a movie of worm with motion back in
    for ii = 1:1:size(Xrecon,1) % only plotting every 10 frames to speed up movie
%   
        subplot(5,4,[13,14])
        imshow(reshape(c_L{1}(:,ii),13,20))
        caxis([-3,3])
        colorbar
        subplot(5,4,[15,16])
        imshow(reshape(c_L{2}(:,ii),13,20))
        caxis([-1.5,1.5])
        colorbar
        subplot(5,4,[17,18])
        imshow(reshape(c_L{3}(:,ii),13,20))
        caxis([-1.5,1.5])
        colorbar
        subplot(5,4,[19,20])
        imshow(reshape(c_L{4}(:,ii),13,20))  
        caxis([-1.5,1.5])
        colormap jet
        colorbar
        
        color_vec = colormap;
        subplot(5,4,[1:3,5:7,9:11])
        plot(Xrecon(ii, :)- Xrecon(1, 1), Yrecon(ii, :) - Yrecon(1, 1), 'Color', color_vec(round(255*((neu_val+3)/6)),:),'LineWidth',3)
        hold on
        plot(Xrecon(ii,1)- Xrecon(1, 1), Yrecon(ii, 1) - Yrecon(1, 1), '.', 'Color', 'b', 'MarkerSize', 20)
        axis equal
        xlim([ min(min(Xrecon - Xrecon(1, 1))), max(max(Xrecon - Xrecon(1, 1)))])
        ylim([ min(min(Yrecon - Yrecon(1, 1))), max(max(Yrecon - Yrecon(1, 1)))])
        hold off
        
        cc1=['Index of the neuron is (',num2str(idx),')'];
        text(460,20,cc1);
        neu_val = c_L{idx(1)}((idx(1)-1)*260+idx(2),ii);
        cc2=['The neuron value is: ',num2str(neu_val)];
        text(460,800,cc2); 
        subplot(5,4,8)
        imshow(c_L{idx(1)}((idx(1)-1)*260+idx(2),ii))
        caxis([-3,3])
        colormap jet
        
      %  pause(0.01) % forces refresh so you can see the plot
    end
    
 %% find correlation to features
 fea_info = h5info([cur_folder,feature_name]); 
 fea_timeseries = h5read([cur_folder,feature_name],'/features_timeseries'); 
 SNames = fieldnames(fea_timeseries);
 
 
 
end