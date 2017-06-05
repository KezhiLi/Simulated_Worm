%
%
% 
% 
% 
% 

% c matrix of different layers
c_L = cell(4,1);
y_test = h5read(hdf5_file_name,'/y_test');
c_L{1} = h5read(hdf5_file_name,'/c_L0_col');
c_L{2} = h5read(hdf5_file_name,'/c_L1_col');
c_L{3} = h5read(hdf5_file_name,'/c_L3_col');
c_L{4} = h5read(hdf5_file_name,'/c_L5_col');
eig_vec = h5read(hdf5_file_name,'/eig_vec');
len_vec = h5read(hdf5_file_name,'/len_vec');
%nan_idx = h5read(hdf5_file_name,'/nan_idx');

%% plot
% choose the layer here
layer = 2; 
% draw the time series
figure, imagesc(c_L{layer}),colorbar

