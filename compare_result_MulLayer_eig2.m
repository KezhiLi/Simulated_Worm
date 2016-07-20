addpath('X:\Kezhi\fastICA');


fps = 5;


csv_predicted = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_predicted_full2.csv';
predicted = csvread(csv_predicted);
csv_y_test = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_y_test_full2.csv';
y_test = csvread(csv_y_test);
csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_MulLayer_generated_full2.csv';
generated_ske = csvread(csv_generated);

load_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
eig_vec_file = 'eig_para_full2.hdf5';
eig_vec = h5read([load_path,eig_vec_file],'/eig_vec'); 
%eig_vec = A{n};
len_vec = h5read([load_path,eig_vec_file],'/len_vec'); 

for tt = 3;
    if tt==1
        eig_radias_vec = predicted(:,1:end-1);
        mean_angle_vec = predicted(:,end);
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_12.gif';
    elseif tt ==2
    %    eig_radias_vec = icasig{n};
        eig_radias_vec = y_test(:,1:end-1);
        mean_angle_vec = y_test(:,end);
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_22.gif';
    elseif tt == 3
        eig_radias_vec = generated_ske(:,1:end-1);
        mean_angle_vec = generated_ske(:,end);
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_32.gif';
    end
   
    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end
    
    radias_vec=eig_vec*eig_radias_vec + (kron(ones(1,size(eig_vec,1)), mean_angle_vec))';

    % radias to ske
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
        %plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
        plot(pred_ske(1,:,ii),pred_ske(2,:,ii),'*-');
        axis equal
        pause(0.2)
       % mov(ii) = save_crt_fra(filename_gif,ii, fps);
    end
end