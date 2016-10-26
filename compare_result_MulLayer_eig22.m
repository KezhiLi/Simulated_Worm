% 6 dimensions + 1 DiffOfMeanofAbsAngles


addpath('X:\Kezhi\fastICA');


fps = 5;


csv_predicted = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_predicted_full22.csv';
predicted = csvread(csv_predicted);
csv_y_test = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_y_test_full22.csv';
y_test = csvread(csv_y_test);
csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_generated_full22.csv';
generated_ske = csvread(csv_generated);

load_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
eig_vec_file = 'eig_para_full22.hdf5';
eig_vec = h5read([load_path,eig_vec_file],'/eig_vec'); 
%eig_vec = A{n};
len_vec = h5read([load_path,eig_vec_file],'/len_vec'); 

FirstNoFrm = length(mean_angle_vec)-size(y_test,1);
FirstAbsAng = mean_angle_vec(FirstNoFrm);

for tt = 3;
    if tt==1
        eig_radias_vec = predicted(:,1:end-1);
        mean_angle_vec_diff = predicted(:,end);
        mean_angle_vec_cur = (mean_angle_vec(FirstNoFrm:end-1))'+mean_angle_vec_diff;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_1(22).gif';
        
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
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_2(22).gif';
    elseif tt == 3
        eig_radias_vec = generated_ske(:,1:end-1);
        mean_angle_vec_diff = generated_ske(:,end);
        mean_angle_vec_diff(1) = mean_angle_vec_diff(1) + FirstAbsAng;
        mean_angle_vec_cur = cumsum(mean_angle_vec_diff);
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_3(22).gif';
    end
   
    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end
    
    radias_vec=eig_vec*eig_radias_vec + (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur))';

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
    
    if tt == 1
        if size(eig_radias_vec2,1)>size(eig_radias_vec2,2)
            eig_radias_vec2 = eig_radias_vec2';
        end

        radias_vec2=eig_vec*eig_radias_vec2 + (kron(ones(1,size(eig_vec,1)), mean_angle_vec_cur2))';

        % radias to ske
        rho2 = median(len_vec,2);
        pred_ske_diff2 = zeros(2,size(radias_vec2,1),size(radias_vec2,2));
        pred_ske2 = zeros(2,size(radias_vec2,1)+1,size(radias_vec2,2));
        for ii = 1:size(radias_vec2,2);
            [pred_ske_diff2(1,:,ii), pred_ske_diff2(2,:,ii)] = pol2cart(radias_vec2(:,ii),rho2);
            pred_ske2(:,2:end,ii) = cumsum(pred_ske_diff2(:,:,ii),2);
            oriPoint2 = floor(size(pred_ske2,2)+1)/2;
            pred_ske2(:,:,ii) = pred_ske2(:,:,ii) - kron(ones(1,size(pred_ske2,2)),pred_ske2(:,oriPoint2,ii));
        end        
        
    end
    
        
    % show the animation of skeleton 
    for ii = 1:size(pred_ske,3);
        %img = zeros(480,640);
        %plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
        plot(pred_ske(1,:,ii),pred_ske(2,:,ii),'r*-');
        if tt ==1
            hold on, 
            % real skeleton
            plot(pred_ske2(1,:,ii),pred_ske2(2,:,ii),'g*-');
            % last skeleton
            if ii>1
                plot(pred_ske(1,:,ii-1),pred_ske2(2,:,ii-1),'y*-');
            end
            hold off,
        end
       axis equal
        pause(0.15)
   %     mov(ii) = save_crt_fra(filename_gif,ii, fps);
    end
end