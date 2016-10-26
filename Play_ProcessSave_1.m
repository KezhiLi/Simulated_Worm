% 5 dimensions + 1 DiffOfMeanofAbsAngles
% play the predicted videos


%%
% read_pcaSke_full21
% 5 eigen coefficients
% 1 difference of mean of the absolute angle (after adjusted)

addpath('X:\Kezhi\fastICA');
cur_path = 'C:\Users\kezhili\Documents\Python Scripts\';
load_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
save_root = 'Z:\DLWeights\';

% find all 'on food' file names
all_file_incSwim = subdir('Z:\Results\nas207-1\*_skeletons.hdf5');
num_all_file = size(all_file_incSwim,1);
NonSwimInd = [];
jj = 1;
for ii = 1:num_all_file;
    if isempty(regexp(all_file_incSwim(ii).name,'swimming'))&&~isempty(regexp(all_file_incSwim(ii).name,'food'))
        NonSwimInd = [NonSwimInd, ii];
    end
end
all_file = all_file_incSwim(NonSwimInd);

% hdf5_file_name = 'Z:\Results\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11_skeletons.hdf5';
% masked_image_file =  'Z:\MaskedVideos\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\800 MY16 on food L_2011_03_03__16_52___3___11.hdf5';
% file_name = '247 JU438 on food L_2011_03_03__11_18___3___1';
% hdf5_file_name = ['Z:\Results\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\',file_name,'_skeletons.hdf5'];

time_step = 3;

num_file = 10; %1096
tt = 0;
window_len =50;

% run all files
for nf = 476:476;  % 476
    disp([num2str(nf),'/',num2str(num_file)])
    hdf5_file_name = all_file(nf).name
    % obtain the filename
    [pathstr,file_name,ext]=fileparts(hdf5_file_name);

%    try
        save_path = [save_root,pathstr(12:end),'\',file_name(1:end-10)];
        
        % read information
        ske_info = h5info(hdf5_file_name, '/skeleton');
        frame_size = ske_info.Dataspace.Size(1:2);
        frame_total = ske_info.Dataspace.Size(3);
        ske_full = h5read(hdf5_file_name,'/skeleton'); 


        % subsample skeleton
        ske = ske_full(:,1:time_step:end,1:time_step:end);
        ske_diff = ske(:,2:end,:)-ske(:,1:end-1,:);

        % handle NaN
        cut_skediff = ske_diff;
        cut_skediff_ind = 1;
        sample_skediff =[];
        for ii = 1:size(cut_skediff ,3)
            if ~isnan(sum(sum(cut_skediff(:,:,ii))))
                sample_skediff(:,:,cut_skediff_ind) = cut_skediff(:,:,ii);
                cut_skediff_ind = cut_skediff_ind+1;
            end
        end
        % initialize angle vectors
        angle_vec = zeros(size(sample_skediff,2),size(sample_skediff,3));
        % initialize length of vectors
        len_vec = angle_vec;
        % cart cooridnates to pollar
        for ii = 1:size(sample_skediff,3);
            [angle_vec(:,ii), len_vec(:,ii)] = cart2pol(sample_skediff(1,:,ii),sample_skediff(2,:,ii));
        end

        angle_vec_ori = angle_vec;

        % make angle_vec continious
        for jj = 1:size(angle_vec,2)-1;
            jj
            while angle_vec(1,jj)- angle_vec(1,jj+1)>pi 
                angle_vec(1,jj+1) = 2*pi + angle_vec(1,jj+1);
            end
            while  angle_vec(1,jj)- angle_vec(1,jj+1) < -pi
                angle_vec(1,jj+1) = -2*pi + angle_vec(1,jj+1);
            end
            for ii = 1:size(angle_vec,1)-1;
                while angle_vec(ii,jj)- angle_vec(ii+1,jj)>pi 
                    angle_vec(ii+1,jj) = 2*pi + angle_vec(ii+1,jj);
                end
                while angle_vec(ii,jj)- angle_vec(ii+1,jj)< - pi 
                    angle_vec(ii+1,jj) = -2*pi + angle_vec(ii+1,jj);
                end
            end
        end
        % cope with the last column
        jj = size(angle_vec,2);
        for ii = 1:size(angle_vec,1)-1;
            while angle_vec(ii,jj)- angle_vec(ii+1,jj)>pi 
                angle_vec(ii+1,jj) = 2*pi + angle_vec(ii+1,jj);
            end
            while angle_vec(ii,jj)- angle_vec(ii+1,jj)< - pi 
                angle_vec(ii+1,jj) = -2*pi + angle_vec(ii+1,jj);
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

        for i=1:n
            [icasig{i,1}, A{i,1}, W{i,1}]=fastica(angle_vec_adj,'approach','defl', ...
                'numOfIC', lasteig(1,i), 'g', 'pow3', 'finetune', 'pow3', ...
                'stabilization', 'on', 'lasteig', lasteig(1,i));
        end

        eig_coef = [icasig{n};diff_mean_angle_vec];
        eig_vec = A{n};

        %%
         radias_vec=A{n}*icasig{n}; % + (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur))';
         %radias_vec_shouldbe= angle_vec;
         %radias_vec_shouldbe= A{n}*icasig{n} + kron(ones(size(angle_vec,1),1),mean_angle_vec);
         %radias_vec_shouldbe=  angle_vec_adj + kron(ones(size(angle_vec,1),1),mean_angle_vec);
         radias_vec_shouldbe=  angle_vec_adj;
         
%          for ii = 1:size(pred_ske,3);
%              ii
%             %img = zeros(480,640);
%             %plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
%             plot(radias_vec(:,ii),'r*-');
%             hold on,
%             plot(radias_vec_shouldbe(:,ii),'g*-');
%             hold off,
%             axis equal
%             pause(0.15)
%           %  mov(ii) = save_crt_fra(filename_gif,ii, fps);
%          end
        %% generate data

        % save_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
        
%         if ~exist([save_root,pathstr(12:end)], 'dir')
%             mkdir([save_root,pathstr(12:end)]);
%         end
% 
%         if exist([save_path,'_eig.hdf5'])==2
%             delete([save_path,'_eig.hdf5']);
%         end
%         h5create([save_path,'_eig.hdf5'], '/eig_coef', size(eig_coef), 'Datatype', 'double', ...
%             'Chunksize', size(eig_coef), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
%         h5write([save_path,'_eig.hdf5'], '/eig_coef', eig_coef);
%         h5create([save_path,'_eig.hdf5'], '/eig_vec', size(eig_vec), 'Datatype', 'double', ...
%             'Chunksize', size(eig_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
%         h5write([save_path,'_eig.hdf5'], '/eig_vec', eig_vec);
%         h5create([save_path,'_eig.hdf5'], '/len_vec', size(len_vec), 'Datatype', 'double', ...
%             'Chunksize', size(len_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
%         h5write([save_path,'_eig.hdf5'], '/len_vec', len_vec);
%         h5create([save_path,'_eig.hdf5'], '/mean_angle_vec', size(mean_angle_vec), 'Datatype', 'double', ...
%             'Chunksize', size(mean_angle_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
%         h5write([save_path,'_eig.hdf5'], '/mean_angle_vec', mean_angle_vec);
%         
        %%

        %load_path = 'Z:\DLWeights\nas207-1\experimentBackup\from pc207-18\!worm_videos\copied from pc207-13\misc_videos\rb1546\';
        %eig_vec_file = 'RB1546 on food L_2010_11_11__12_44_04___1___6_eig';
        %csv_generated = [save_path,'_y_test-hid2-2.csv'];
        %csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_generated_full21(247JU438).csv';
        csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_generated_full22.csv';
        
        %csv_generated = [save_path,'_eig(4hid).csv']; %
        %csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\generated.csv';
        %csv_generated = 'Z:\DLWeights_temp\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\03-03-11\247 JU438 on food L_2011_03_03__11_18___3___1_generated-hid3-4.csv';
        
        generated_ske = csvread(csv_generated);
     %   generated_ske = eig_coef(:,(round(0.9*size(eig_coef,2))):(round(0.9*size(eig_coef,2)))+400)';

%          eig_vec = h5read([save_path,'_eig.hdf5'],'/eig_vec');
%          len_vec = h5read([save_path,'_eig.hdf5'],'/len_vec');
%          mean_angle_vec = h5read([save_path,'_eig.hdf5'],'/mean_angle_vec');

        FirstNoFrm = round(length(mean_angle_vec)*0.9+window_len);
        FirstAbsAng = mean_angle_vec(FirstNoFrm);

        eig_radias_vec = generated_ske(:,1:end-1);
        mean_angle_vec_diff = generated_ske(:,end);
        mean_angle_vec_diff(1) = mean_angle_vec_diff(1) + FirstAbsAng;
        mean_angle_vec_cur = cumsum(mean_angle_vec_diff);
        filename_gif = [save_path,'_eig.gif'];


        if size(eig_radias_vec,1)>size(eig_radias_vec,2)
            eig_radias_vec = eig_radias_vec';
        end

          radias_vec=eig_vec*eig_radias_vec + (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur))';
%          %radias_vec_shouldbe= angle_vec;
%          %radias_vec_shouldbe= A{n}*icasig{n} + kron(ones(size(angle_vec,1),1),mean_angle_vec);
%          radias_vec_shouldbe=  angle_vec(:,round(0.9*size(angle_vec,2)+window_len):(round(0.9*size(angle_vec,2))+window_len+400)); % + kron(ones(size(angle_vec,1),1),mean_angle_vec);

         
%          for ii = 1:size(pred_ske,3);
%             %img = zeros(480,640);
%             %plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
%             plot(radias_vec(:,ii),'r*-');
%             hold on,
%             plot(radias_vec_shouldbe(:,ii),'g*-');
%             hold off,
%             axis equal
%             pause(0.15)
%           %  mov(ii) = save_crt_fra(filename_gif,ii, fps);
%           end
        
                 
         
        % radias to ske
        rho = median(len_vec,2);
        pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
        pred_ske = zeros(2,size(radias_vec,1)+1,size(radias_vec,2));
        for ii = 1:size(radias_vec,2);
            [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho);
            pred_ske(:,2:end,ii) = cumsum(pred_ske_diff(:,:,ii),2);
            oriPoint = floor(size(pred_ske,2)+1)/2;
            pred_ske(:,:,ii) = pred_ske(:,:,ii) - kron(ones(1,size(pred_ske,2)),pred_ske(:,oriPoint,ii));
        end

%         rho = median(len_vec,2);
%         pred_ske_diff_shouldbe = zeros(2,size(radias_vec_shouldbe,1),size(radias_vec_shouldbe,2));
%         pred_ske_shouldbe = zeros(2,size(radias_vec_shouldbe,1)+1,size(radias_vec_shouldbe,2));
%         for ii = 1:size(radias_vec_shouldbe,2);
%             [pred_ske_diff_shouldbe(1,:,ii), pred_ske_diff_shouldbe(2,:,ii)] = pol2cart(radias_vec_shouldbe(:,ii),rho);
%             pred_ske_shouldbe(:,2:end,ii) = cumsum(pred_ske_diff_shouldbe(:,:,ii),2);
%             oriPoint = floor(size(pred_ske_shouldbe,2)+1)/2;
%             pred_ske_shouldbe(:,:,ii) = pred_ske_shouldbe(:,:,ii) - kron(ones(1,size(pred_ske_shouldbe,2)),pred_ske_shouldbe(:,oriPoint,ii));
%         end
%         
        % show the animation of skeleton
        for ii = 1:size(pred_ske,3);
            %img = zeros(480,640);
            %plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
            plot(pred_ske(1,:,ii),pred_ske(2,:,ii),'r*-');
%             hold on,
%             plot(pred_ske_shouldbe(1,:,ii),pred_ske_shouldbe(2,:,ii),'g*-');
%             hold off,
            axis equal
            pause(0.15)
          %  mov(ii) = save_crt_fra(filename_gif,ii, fps);
        end


%    catch ME
        tt = tt+1;
%    end
end
