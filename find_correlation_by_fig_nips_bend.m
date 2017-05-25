% find correlation by figures
windowSize = 150;
b = (1/windowSize)*ones(1,windowSize);
a = 1;


startpt = 1300;% 1 
endpt =  2500;% length(comb_mtx_NaN0{layer}(feature,:)); 

layer= 4;
feature = 13; % 39: middle speed; 13: midbody bend
for ii =  71 % 57:1:316;
    feature_vec = comb_mtx_NaN0{layer}(feature,:);
    feature_vec = feature_vec/norm(feature_vec);
    feature_vec = feature_vec-mean(feature_vec);
    figure, plot(feature_vec(startpt:endpt),'r'),   % 3050:6700
    hold on,   
%    plot(comb_mtx_NaN0{1}(43,5000:6000),'m'), 
    neuron_vec = comb_mtx_NaN0{layer}(ii,:);
    neuron_vec = neuron_vec/norm(neuron_vec);
    neuron_vec = neuron_vec - mean(neuron_vec);
    plot(neuron_vec(startpt:endpt));
    plot(zeros(1,length(neuron_vec(startpt:endpt))),'k--')
    %plot(neuron_vec*sqrt(norm(comb_mtx_NaN0{layer}(feature,:)))) %3050:6700
    axis([0 1200  -0.024 0.024])
   
    title('bend vs. memory')
    xlabel('frame indexes')
    ylabel('normalized magnitudes')
    legend('bend','memory')

%     y0 = filter(b,a,feature_vec); %3050:6700 
%     y = filter(b,a,neuron_vec);%3050:6700
% 
%     figure, plot(y0(startpt:endpt),'r'), 
%     hold on,  
%     plot(y(startpt:endpt))
%     plot(zeros(1,length(y(startpt:endpt))),'k--')
%     axis([0 1200  -0.0125 0.01  ])
end

% figure, subplot(2,1,1)
% hold on
% plot(feature_vec(startpt:endpt),'r'),
% plot(neuron_vec(startpt:endpt));
% plot(zeros(1,length(neuron_vec(startpt:endpt))),'k--')
% axis([0 1200  -0.05 0.022  ])
% title('speed vs. memory')
% xlabel('frame indexes')
% ylabel('normalized magnitudes')
% legend('speed','memory')
% subplot(2,1,2)
% plot(y0(startpt:endpt),'r'), 
% hold on,
% plot(y(startpt:endpt))
% plot(zeros(1,length(y(startpt:endpt))),'k--')
% axis([0 1200  -0.014 0.012  ])
% title('filtered speed vs. filtered memory')
% xlabel('frame indexes')
% ylabel('normalized magnitudes')
% legend('filtered speed','filtered memory')
% 
% cor1 = corr(feature_vec', neuron_vec');
% cor2 = corr(y0', y');
