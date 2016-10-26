% find all .fig file names
all_fig_file = subdir('C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\7-200-200-200-200-7\07-03-11\600epochs\no_noise\');
%'C:\Users\kezhili\Documents\Python Scripts\data\FromAWS\08-03-11\1000epochs\*.fig'

num_file = size(all_fig_file,1);

 
for nf = 1:num_file;  % 476
    disp([num2str(nf),'/',num2str(num_file)])
    fig_file_name = all_fig_file(nf).name
    open(fig_file_name);
    print('-dpng', '-r600', strcat(fig_file_name(1:end-4),'.png'));
    close all
end