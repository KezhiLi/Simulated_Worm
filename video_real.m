% function: produce sample real video clips for Turing test
%
%
%
%
%  15/03/2017 Kenneth (Kezhi) Li @Imperial College, CSC, MRC


addpath('X:\Kezhi\fastICA');
addpath('C:\Kezhi\MyCode!!!\Simulated_Worm\Test_Eric\');

input_folder = ['Z:\DLWeights\nas207-1\experimentBackup\*_eig.hdf5'];
all_file = subdir(input_folder);
num_file = length(all_file);  


for nf = 1013:num_file  %num_csv;  % 476

    disp([num2str(nf),'/',num2str(num_file)])
    hdf5_file_name = all_file(nf).name
    % obtain the filename
    [pathstr,file_name,ext]=fileparts(hdf5_file_name);
    fname_head = ['Z:\Ken_Samples\real\',pathstr(14:end),'\',file_name(1:end-3)];
    
%     if exist(['Z:\Ken_Samples\real\',pathstr(14:end),'\'])
%         continue;
%     end
    
    % read information
    eig_coef = h5read(hdf5_file_name,'/eig_coef');
    eig_vec = h5read(hdf5_file_name,'/eig_vec');    
    len_vec = h5read(hdf5_file_name,'/len_vec');  
    frame_total = size(eig_coef,2);
    
    eig_radias_vec = eig_coef(1:size(eig_vec,2),:);
    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end
    
    radias_vec=eig_vec*eig_radias_vec; 
    
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
    frm_len = 600;
    if frame_total<frm_len
        start_frm = 1;
        end_frm = frame_total;
    else
        start_frm = max(1,round(rand(1)*(frame_total- frm_len-10)));
        end_frm = min((start_frm+frm_len),frame_total);
    end
    
    fname = [fname_head,'samp_',num2str(start_frm),'.avi'];
    [pathstr2,file_name2,ext2]=fileparts(fname);
    if ~exist(pathstr2)
        mkdir(pathstr2);
    elseif exist(fname)
        delete(fname)
    end
    aviobj = VideoWriter(fname);
    aviobj.FrameRate = fps;
    open(aviobj)
    
    % cut Xrecon,Yrecon
    Xrecon = Xrecon(start_frm: end_frm,:);
    Yrecon = Yrecon(start_frm: end_frm,:);
    % generate video
    for ii = 1: size(Xrecon,1) %1: size(Xrecon,1) % only plotting every 10 frames to speed up movie 
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
end
