% function: show intermediam c matrices for a worm video
%
%
%
%
%  14/02/2017 Kezhi Li @Imperial College, CSC, MRC

addpath('X:\Kezhi\fastICA');
addpath('C:\Kezhi\MyCode!!!\Simulated_Worm\Test_Eric\');

% strain name
stain =  'N2';
% the hdf5 file storages the c matrices
%strain_name =  'N2 on food L_2011_03_24__12_06_06___1___1';
strain_name = 'N2 on food L_2011_03_17__11_34_06___7___3';  % this one good

% current folder
cur_folder = ['Z:\DLWeights\eig_catagory_Straits\',stain,'\',strain_name,'\'];
%'N2 on food L_2011_03_17__11_34_06___7___3';
name = ['interm_',strain_name,'_eig_7-260-260-260-260-7_600ep.hdf5'];
feature_name = [strain_name,'_features.hdf5'];
hdf5_file_name = [cur_folder,name];


idx_layer = 1;
idx_neuron = 120 ;
idx_feature = 39;

% the index of layer and neuron one wants to observe
idx = [idx_layer, idx_neuron];

% c matrix of different layers
c_L = cell(4,1);
res_L = cell(4,1);
y_test = h5read(hdf5_file_name,'/y_test');
c_L{1} = h5read(hdf5_file_name,'/c_L0_col');
c_L{2} = h5read(hdf5_file_name,'/c_L1_col');
c_L{3} = h5read(hdf5_file_name,'/c_L3_col');
c_L{4} = h5read(hdf5_file_name,'/c_L5_col');
res_L{1} = h5read(hdf5_file_name,'/res_L0_col');
res_L{2} = h5read(hdf5_file_name,'/res_L1_col');
res_L{3} = h5read(hdf5_file_name,'/res_L3_col');
res_L{4} = h5read(hdf5_file_name,'/res_L5_col');
eig_vec = h5read(hdf5_file_name,'/eig_vec');
len_vec = h5read(hdf5_file_name,'/len_vec');
nan_idx = h5read(hdf5_file_name,'/nan_idx');

%c_L = [c_L0_col;c_L1_col;c_L3_col;c_L5_col];




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
    
%     % plot a movie of worm with motion back in
%     for ii = 3000:1:3800  %1: size(Xrecon,1) % only plotting every 10 frames to speed up movie 
%         
%                 subplot(5,4,[13,14])
%                 imshow(reshape(c_L{1}(:,ii),13,20))
%                 title('1st layer')
%                 caxis([-3,3])
%                 colorbar
%                 subplot(5,4,[15,16])
%                 imshow(reshape(c_L{2}(:,ii),13,20))
%                 title('2nd layer')
%                 caxis([-1.5,1.5])
%                 colorbar
%                 subplot(5,4,[17,18])
%                 imshow(reshape(c_L{3}(:,ii),13,20))
%                 title('3rd layer')
%                 caxis([-1.5,1.5])
%                 colorbar
%                 subplot(5,4,[19,20])
%                 imshow(reshape(c_L{4}(:,ii),13,20))
%                 title('4th layer')
%                 caxis([-1.5,1.5])
%         colormap jet
%                 colorbar
%         
% %                 hp4 = get(subplot(5,4,20),'Position');
% %                 colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*1.1])
%         
% 
%         
%         color_vec = colormap;
%         subplot(5,4,[1:3,5:7,9:11])
%         neu_val = min(3,max(-2.98,c_L{idx(1)}((idx(1)-1)*260+idx(2),ii)));
%         plot(Xrecon(ii, :)- Xrecon(1, 1), Yrecon(ii, :) - Yrecon(1, 1), 'Color', color_vec(round(255*((neu_val+3)/6)),:),'LineWidth',3)
%         hold on
%         plot(Xrecon(ii,1)- Xrecon(1, 1), Yrecon(ii, 1) - Yrecon(1, 1), '.', 'Color', 'b', 'MarkerSize', 20)
%         axis equal
%         xlim([ min(min(Xrecon - Xrecon(1, 1))), max(max(Xrecon - Xrecon(1, 1)))])
%         ylim([ min(min(Yrecon - Yrecon(1, 1))), max(max(Yrecon - Yrecon(1, 1)))])
%         hold off
%         
%         fea_val = fea_timeseries_show(idx_feature,ii);
%         cc0 = ['The feature''s value is: ',num2str(fea_val)];
%         text(1860,1250,cc0);
%         cc1=['Index of the neuron is (',num2str(idx),')'];
%         text(1860,-400,cc1);
%         cc2=['The neuron''c value is: ',num2str(neu_val)];
%         text(1860,1100,cc2);
%         subplot(5,4,8)
%         imshow(c_L{idx(1)}((idx(1)-1)*260+idx(2),ii))
%         caxis([-3,3])
%         colormap jet
%         
%         pause(0.00000000001) % forces refresh so you can see the plot
%     end
    
end

%% find correlation to features
fea_info = h5info([cur_folder,feature_name]);
fea_timeseries = h5read([cur_folder,feature_name],'/features_timeseries');
SNames = fieldnames(fea_timeseries);

%  midbody_bends_length = length(worm.locomotion.bends.midbody.amplitude);
%  midbody_bends = worm.locomotion.bends.midbody.amplitude(round(0.1*midbody_bends_length):3:end);

fea_timeseries_all = zeros(length(SNames),max(size(fea_timeseries.(SNames{1}))));
for ii = 1:length(SNames);
    fea_timeseries_all(ii,:) = fea_timeseries.(SNames{ii});
end
subsample_step = 3;
fea_timeseries_subs = fea_timeseries_all(:,1:subsample_step:end); 
% abandon the NaN that have been abandoned previously
fea_timeseries_subs_noNaN = fea_timeseries_subs(:,setdiff(1:size(fea_timeseries_subs,2),nan_idx));
% find the start index
idx_timeseries = size(fea_timeseries_subs_noNaN,2)- size(c_L{1},2)+1;
%idx_timeseries = length(fea_timeseries_subs_noNaN)*0.1+50;
fea_timeseries_show =  fea_timeseries_subs_noNaN(:,idx_timeseries:end); 
%fea_timeseries_show = fea_timeseries_subs()
row_num = size(fea_timeseries_show,1);
lay_num = size(c_L,1);
idx_num = size(c_L{1},1);
col_num = lay_num*idx_num;
corr_mtx = zeros(row_num,col_num);
%  for ii = 1:row_num
%      for jj = 1:lay_num
%          for kk = 1:idx_num
%             corr_mtx(ii,(jj-1)*idx_num+kk) =xcorr( fea_timeseries_show(ii,:),c_L{jj,1}(kk,:));
%          end
%      end
%  end
comb_mtx = cell(4,1);
comb_mtx_NaN0 = cell(4,1);
R = cell(4,1);
R_largeEntry = cell(4,1);
thre_coe = 0.2;
for ii = 1:4;
    comb_mtx{ii} = [fea_timeseries_show;c_L{ii,1}];
    comb_mtx_NaN0{ii} = comb_mtx{ii};
    comb_mtx_NaN0{ii}(isnan(comb_mtx_NaN0{ii}))=0;
    
    R{ii} = corrcoef(comb_mtx_NaN0{ii}');
    R_largeEntry{ii} = zeros(size(R{ii}));
    R_largeEntry{ii}(abs(R{ii})>thre_coe ) = R{ii}(abs(R{ii})>thre_coe );
    %  R_largeEntry{ii}(1:56,1:56)=0;
    %  R_largeEntry{ii}(57:end,57:end)=0;
    
%     for jj = 1:size(comb_mtx_NaN0{ii},1)
%         comb_mtx_NaN0{ii}(jj,:) = comb_mtx_NaN0{ii}(jj,:)-mean(comb_mtx_NaN0{ii}(jj,:));
%         comb_mtx_NaN0{ii}(jj,:) = comb_mtx_NaN0{ii}(jj,:)/std(comb_mtx_NaN0{ii}(jj,:));
%     end
end

comb_mtx2 = cell(4,1);
comb_mtx2_NaN0 = cell(4,1);
R2 = cell(4,1);
R2_largeEntry = cell(4,1);
for ii = 1:4;
    comb_mtx2{ii} = [fea_timeseries_show;res_L{ii,1}];
    comb_mtx2_NaN0{ii} = comb_mtx2{ii};
    comb_mtx2_NaN0{ii}(isnan(comb_mtx2_NaN0{ii}))=0;
    
    R2{ii} = corrcoef(comb_mtx2_NaN0{ii}');
    R2_largeEntry{ii} = zeros(size(R2{ii}));
    R2_largeEntry{ii}(abs(R2{ii})>thre_coe ) = R2{ii}(abs(R2{ii})>thre_coe );
    %  R_largeEntry{ii}(1:56,1:56)=0;
    %  R_largeEntry{ii}(57:end,57:end)=0;
    
%     for jj = 1:size(comb_mtx2_NaN0{ii},1)
%         comb_mtx2_NaN0{ii}(jj,:) = comb_mtx2_NaN0{ii}(jj,:)-mean(comb_mtx2_NaN0{ii}(jj,:));
%         comb_mtx2_NaN0{ii}(jj,:) = comb_mtx2_NaN0{ii}(jj,:)/std(comb_mtx2_NaN0{ii}(jj,:));
%     end
end

% for ii =1:56;
%     figure, plot(comb_mtx_NaN0{1,1}(57,:),'r'), hold on, plot(comb_mtx_NaN0{1,1}(ii,:))
% end

comb_mtx3 = cell(4,1);
comb_mtx3_NaN0 = cell(4,1);
R3 = cell(4,1);
R3_largeEntry = cell(4,1);
for ii = 1:4;
    comb_mtx3{ii} = [fea_timeseries_show(:,2:end);1e2*abs(res_L{ii}(:,2:end)-res_L{ii}(:,1:end-1))];
    comb_mtx3_NaN0{ii} = comb_mtx3{ii};
    comb_mtx3_NaN0{ii}(isnan(comb_mtx3_NaN0{ii}))=0;
    
    R3{ii} = corrcoef(comb_mtx3_NaN0{ii}');
    R3_largeEntry{ii} = zeros(size(R3{ii}));
    R3_largeEntry{ii}(abs(R3{ii})>thre_coe ) = R3{ii}(abs(R3{ii})>thre_coe );
    %
    %     for jj = 1:size(comb_mtx3_NaN0{ii},1)
    %         comb_mtx3_NaN0{ii}(jj,:) = comb_mtx3_NaN0{ii}(jj,:)-mean(comb_mtx3_NaN0{ii}(jj,:));
    %         comb_mtx3_NaN0{ii}(jj,:) = comb_mtx3_NaN0{ii}(jj,:)/std(comb_mtx3_NaN0{ii}(jj,:));
    %     end
end