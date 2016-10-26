
fps = 5;

csv_predicted2 = 'C:\Users\kezhili\Documents\Python Scripts\data\predicted2_full(18-500-500-18).csv';
predicted2 = csvread(csv_predicted2);

predicted = predicted2(:,1:end-2);
predicted_loc_diff = predicted2(:,end-1:end);
predicted_loc = cumsum(predicted_loc_diff);
predicted_loc(1,:) = 0; 


csv_y_test2 = 'C:\Users\kezhili\Documents\Python Scripts\data\y_test2_full(18-500-500-18).csv';
y_test2 = csvread(csv_y_test2);

y_test = y_test2(:,1:end-2);
y_test_loc_diff = y_test2(:,end-1:end);
y_test_loc = cumsum(y_test_loc_diff);
y_test_loc(1,:) = 0; 


csv_generated2 = 'C:\Users\kezhili\Documents\Python Scripts\data\generated2_full(18-500-500-18).csv';
generated_ske2 = csvread(csv_generated2);

generated_ske = generated_ske2(:,1:end-2);
generated_loc_diff = generated_ske2(:,end-1:end);
generated_loc = cumsum(generated_loc_diff);
generated_loc(1,:) = 0; 

for tt = 3;
    if tt==1
        radias_vec = predicted;
        loc = predicted_loc;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\12_full(18-500-500-18).gif';
    elseif tt ==2
        radias_vec = y_test;
        loc = y_test_loc;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\22_full(18-500-500-18).gif';
    elseif tt == 3
        radias_vec = generated_ske;
        loc = generated_loc;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\32_full(18-500-500-18).gif';
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
        ii
        %img = zeros(480,640);
        plot(pred_ske(1,:,ii),pred_ske(2,:,ii),'*-');
        axis equal
      %  axis([2500 4200 5800 7200])
        pause(0.1)
        mov(ii) = save_crt_fra(filename_gif,ii, fps);
    end
end