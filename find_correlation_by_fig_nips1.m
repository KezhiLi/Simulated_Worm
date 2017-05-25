% find correlation by figures
windowSize = 150;
b = (1/windowSize)*ones(1,windowSize);
a = 1;


startpt = 5250; %3400
endpt = 6450; %4600

layer= 1;
feature = 39; % 39: middle speed; 13: midbody bend
for ii =  269 % 57:1:316;
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
    axis([0 1200  -0.05 0.021  ])
   
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

    y0 = filter(b,a,feature_vec); %3050:6700 
    y = filter(b,a,neuron_vec);%3050:6700

    figure, plot(y0(startpt:endpt),'r'), 
    hold on,  
    plot(y(startpt:endpt))
    plot(zeros(1,length(y(startpt:endpt))),'k--')
    axis([0 1200  -0.0125 0.01  ])
end

% figure, plot(R_largeEntry{1}(13,:),'r')
% hold on, plot(R_largeEntry{4}(13,:),'b')

figure, subplot(2,1,1)
hold on
plot(feature_vec(startpt:endpt),'r'),
plot(neuron_vec(startpt:endpt));
plot(zeros(1,length(neuron_vec(startpt:endpt))),'k--')
axis([0 1200  -0.05 0.022  ])
title('speed vs. memory')
xlabel('frame indexes')
ylabel('normalized magnitudes')
legend('speed','memory')
subplot(2,1,2)
plot(y0(startpt:endpt),'r'), 
hold on,
plot(y(startpt:endpt))
plot(zeros(1,length(y(startpt:endpt))),'k--')
axis([0 1200  -0.014 0.012  ])
title('filtered speed vs. filtered memory')
xlabel('frame indexes')
ylabel('normalized magnitudes')
legend('filtered speed','filtered memory')

cor1 = corr(feature_vec', neuron_vec');
cor2 = corr(y0', y');
