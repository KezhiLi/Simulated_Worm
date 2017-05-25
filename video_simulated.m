% function: produce sample simulated video clips for Turing test
%
%
%
%
%  27/03/2017 Kenneth (Kezhi) Li @Imperial College, CSC, MRC


addpath('X:\Kezhi\fastICA');
addpath('C:\Kezhi\MyCode!!!\Simulated_Worm\Test_Eric\');

input_folder = 'Z:\DLWeights\nas207-1\experimentBackup\*_eig.hdf5';
all_file = subdir(input_folder);
num_file = length(all_file);  

%% loop to generate the frames
for nf = 1:num_file  %num_csv;  % 476

    disp([num2str(nf),'/',num2str(num_file)])
    hdf5_file_name = all_file(nf).name
    % obtain the filename
    [pathstr,file_name,ext]=fileparts(hdf5_file_name);
    try
    if ~(~isempty(strfind(file_name,'N2'))||~isempty(strfind(file_name,'unc-8'))||~isempty(strfind(file_name,'ser-6'))||...
            ~isempty(strfind(file_name,'tdc-1'))||~isempty(strfind(file_name,'tbh-1'))||~isempty(strfind(file_name,'cb4856'))||...
            ~isempty(strfind(file_name,'unc-9'))||~isempty(strfind(file_name,'trp-4'))||~isempty(strfind(file_name,'MY'))||...
        ~isempty(strfind(file_name,'LS'))||~isempty(strfind(file_name,'JU'))||~isempty(strfind(file_name,'ED'))||...
        ~isempty(strfind(file_name,'CB'))||~isempty(strfind(file_name,'AQ')))
        continue;
    end    
    fname_head = ['Z:\Ken_Samples\simulated\',pathstr(14:end),'\',file_name(1:end-3)];
    
%     if exist(['Z:\Ken_Samples\real\',pathstr(14:end),'\'])
%         continue;
%     end

    % read simulated csv files
    csv_generated_samp = ['Z:\Ken_Samples\simulated_csv\',pathstr(14:end),'\',file_name(1:end-4),'_samp_*.csv'];
    csv_generated = subdir(csv_generated_samp);
    [pathstr_csv,file_name_csv,ext_csv]=fileparts(csv_generated(1).name);
    text_idx = file_name_csv((strfind(file_name_csv,'samp_')+5):end);
    generated_ske = csvread(csv_generated(1).name);
    % eigen vectors
    eig_radias_vec = generated_ske(:,1:end-1);
    % set the frame length <= 600
    frm_len = 600;
    if size(eig_radias_vec,1)>frm_len
        eig_radias_vec = eig_radias_vec(1:frm_len,:);
    end
    % read information
    eig_vec = h5read(hdf5_file_name,'/eig_vec');    
    len_vec = h5read(hdf5_file_name,'/len_vec');  
    
    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end
    
    frame_total = size(eig_radias_vec,2);
    radias_vec = eig_vec*eig_radias_vec; 
    radias_vec = radias_vec + ones(size(radias_vec))*rand(1)*pi;
    
    % radias to ske
    rho = median(len_vec,2);
    pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
    pred_ske = zeros(2,size(radias_vec,1)+1,size(radias_vec,2));
    for ii = 1:size(radias_vec,2);
        [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho); % len_vec(:,ii)
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
    fps = 10;
    
    fname = [fname_head,'simu_',text_idx,'.avi'];
    [pathstr2,file_name2,ext2]=fileparts(fname);
    if ~exist(pathstr2)
        mkdir(pathstr2);
    elseif exist(fname)
        delete(fname)
    end
    aviobj = VideoWriter(fname);
    aviobj.FrameRate = fps;
    open(aviobj)
    
    for ii = 1: frame_total %1: size(Xrecon,1) % only plotting every 10 frames to speed up movie 
        fig = figure;
        plot(Xrecon(ii, :)- Xrecon(1, 1), Yrecon(ii, :) - Yrecon(1, 1), 'Color', 'r','LineWidth',3)
        hold on
        plot(Xrecon(ii,1)- Xrecon(1, 1), Yrecon(ii, 1) - Yrecon(1, 1), '.', 'Color', 'b', 'MarkerSize', 20)
        axis equal
        xlim([ min(min(Xrecon - Xrecon(1, 1)))-500, max(max(Xrecon - Xrecon(1, 1)))]+500)
        ylim([ min(min(Yrecon - Yrecon(1, 1)))-500, max(max(Yrecon - Yrecon(1, 1)))]+500)
        hold off
        FF = getframe(fig);
        writeVideo(aviobj,FF);
        close(fig);
        %pause(0.1) % forces refresh so you can see the plot
        %mov = save_crt_fra(filename, (ii-start_frm+1), fps);
    end
    close(aviobj);
    %movie2avi(mov, fname, 'compression', 'None', 'fps', fps);
    catch
        continue;
    end
end
