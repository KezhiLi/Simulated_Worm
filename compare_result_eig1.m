
fps = 5;

csv_predicted = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_predicted1.csv';
predicted = csvread(csv_predicted);
csv_y_test = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_y_test1.csv';
y_test = csvread(csv_y_test);
csv_generated = 'C:\Users\kezhili\Documents\Python Scripts\data\eig_generated1.csv';
generated_ske = csvread(csv_generated);

load_path = 'C:\Users\kezhili\Documents\Python Scripts\data\';
eig_vec_file = 'eig_para1.hdf5';
eig_vec = h5read([load_path,eig_vec_file],'/eig_vec'); 
len_vec = h5read([load_path,eig_vec_file],'/len_vec'); 

for tt = 2;
    if tt==1
        eig_radias_vec = predicted;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_1.gif';
    elseif tt ==2
        eig_radias_vec = y_test;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_2.gif';
    elseif tt == 3
        eig_radias_vec = generated_ske;
        filename_gif = 'C:\Kezhi\MyCode!!!\Simulated_Worm\eig_3.gif';
    end
   
    if size(eig_radias_vec,1)>size(eig_radias_vec,2)
        eig_radias_vec = eig_radias_vec';
    end

    radias_vec=eig_vec*eig_radias_vec;

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
        plot(pred_ske(1,:,ii),480-pred_ske(2,:,ii),'*-');
        axis equal
        pause(0.2)
       % mov(ii) = save_crt_fra(filename_gif,ii, fps);
    end
end