
fps = 5;

csv_predicted3 = 'C:\Users\kezhili\Documents\Python Scripts\data\predicted3(18-600-600-18).csv';
predicted3 = csvread(csv_predicted3);

%predicted_diff = [medfilt1(predicted3(:,end),3,50),predicted3(:,1:end-3)];
predicted_diff = [predicted3(:,end),predicted3(:,1:end-3)];
predicted =  cumsum(predicted_diff,2);
predicted_loc_diff = predicted3(:,end-2:end-1);
predicted_loc = cumsum(predicted_loc_diff);
predicted_loc(1,:) = 0; 


csv_y_test3 = 'C:\Users\kezhili\Documents\Python Scripts\data\y_test3(18-600-600-18).csv';
y_test3 = csvread(csv_y_test3);

y_test_diff = [y_test3(:,end),y_test3(:,1:end-3)];
y_test = cumsum(y_test_diff,2);
y_test_loc_diff = y_test3(:,end-2:end-1);
y_test_loc = cumsum(y_test_loc_diff);
y_test_loc(1,:) = 0; 


csv_generated3 = 'C:\Users\kezhili\Documents\Python Scripts\data\generated3(18-600-600-18).csv';
generated_ske3 = csvread(csv_generated3);

generated_ske_diff = [generated_ske3(:,end),generated_ske3(:,1:end-3)];
generated_ske = cumsum(generated_ske_diff,2);
generated_loc_diff = generated_ske3(:,end-2:end-1);
generated_loc = cumsum(generated_loc_diff);
generated_loc(1,:) = 0; 

for tt = 1;
    if tt==1
        radias_vec = predicted;
        loc = predicted_loc;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\13_ful1.gif';
    elseif tt ==2
        radias_vec = y_test;
        loc = y_test_loc;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\23_ful1.gif';
    elseif tt == 3
        radias_vec = generated_ske;
        loc = generated_loc;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\33_ful1.gif';
    end
    
    if size(radias_vec,1)>size(radias_vec,2)
        radias_vec = radias_vec';
    end   

    if size(len_vec,1)~=size(radias_vec,1)
        error('the skeleton vector is not the same')
    end
    
    num_pt_ske = size(radias_vec,1)+1;
    
    mid_ske_ind = round(( num_pt_ske+1)/2);
    % estimate vector length
    rho = median(median(len_vec(~isnan(len_vec)),2));
    % predicted skeleton
    pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
    pred_ske = zeros(2, num_pt_ske,size(radias_vec,2));
    
    
    for ii = 1:size(radias_vec,2);
        [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho);
        pred_ske(:,2:end,ii) = cumsum(pred_ske_diff(:,:,ii),2);
        
%         diff_loc = loc(ii,:)' - pred_ske(:,mid_ske_ind,ii)+sample_ske(:,mid_ske_ind,1);
%         pred_ske(:,:,ii) = pred_ske(:,:,ii) + kron(ones(1, num_pt_ske),diff_loc);
    end
    % show the animation of skeleton 
    for ii = 1: size(pred_ske,3);
        %img = zeros(480,640);
        plot(pred_ske(1,:,ii),pred_ske(2,:,ii),'*-');
        axis equal
   %     axis([2500 4500 5800 7500])
        pause(0.1)
  %      mov(ii) = save_crt_fra(filename_gif,ii, fps);
    end
end