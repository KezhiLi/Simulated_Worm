% find correlation by figures
windowSize = 200;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

layer= 1;
feature = 38; % 39: middle speed; 13: midbody bend
for ii =  4% 57:1:316;
    figure, plot(comb_mtx_NaN0{layer}(feature,:),'r'),   % 3050:6700
    hold on,   
%    plot(comb_mtx_NaN0{1}(43,5000:6000),'m'), 
    plot(comb_mtx_NaN0{layer}(ii,:)) %3050:6700
  %  axis([0,650,-5,5])
   
%     figure, plot(comb_mtx_NaN0{layer}(feature,:),'r'), %4050:4700
%     hold on, 
% %    plot(comb_mtx_NaN0{1}(43,5000:6000),'m'), 
%     plot(comb_mtx_NaN0{layer}(ii+1,:))
%  %   axis([0,650,-5,5])
%     
%     figure, plot(comb_mtx_NaN0{layer}(feature,:),'r'), 
%     hold on, 
% %    plot(comb_mtx_NaN0{1}(43,5000:6000),'m'), 
%     plot(comb_mtx_NaN0{layer}(ii+2,:))
%  %   axis([0,650,-5,5])

    y0 = filter(b,a,comb_mtx_NaN0{layer}(feature,:)); %3050:6700 
    y = filter(b,a,comb_mtx_NaN0{layer}(ii,:));%3050:6700

    figure, plot(y0,'r'), 
    hold on,  
    plot(y)
    plot(zeros(1,length(y)))
end

% figure, plot(R_largeEntry{1}(13,:),'r')
% hold on, plot(R_largeEntry{4}(13,:),'b')
