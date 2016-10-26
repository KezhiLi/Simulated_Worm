
fps = 5;

csv_predicted = 'C:\Users\kezhili\Documents\Python Scripts\data\predicted.csv';
predicted = csvread(csv_predicted);
csv_y_test = 'C:\Users\kezhili\Documents\Python Scripts\data\y_test.csv';
y_test = csvread(csv_y_test);
csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\generated.csv';
generated_ske = csvread(csv_generated);

for tt = 3;
    if tt==1
        radias_vec = predicted;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\1.gif';
    elseif tt ==2
        radias_vec = y_test;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\2.gif';
    elseif tt == 3
        radias_vec = generated_ske;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\3.gif';
    end
    
    if size(radias_vec,1)>size(radias_vec,2)
        radias_vec = radias_vec';
    end   

    if size(len_vec,1)~=size(radias_vec,1)
        error('the skeleton vector is not the same')
    end

    rho = median(len_vec,2);
    pred_ske_diff = zeros(2,size(radias_vec,1),size(radias_vec,2));
    pred_ske = zeros(2,size(radias_vec,1)+1,size(radias_vec,2));
    for ii = 1:size(radias_vec,2);
        [pred_ske_diff(1,:,ii), pred_ske_diff(2,:,ii)] = pol2cart(radias_vec(:,ii),rho);
        pred_ske(:,2:end,ii) = cumsum(pred_ske_diff(:,:,ii),2);
    end
    % show the animation of skeleton 
    for ii = 1:size(pred_ske,3);
        %img = zeros(480,640);
        plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
        axis equal
        %pause(0.1)
        mov(ii) = save_crt_fra(filename_gif,ii, fps);
    end
end