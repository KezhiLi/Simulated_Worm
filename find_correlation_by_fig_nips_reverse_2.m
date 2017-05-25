% find correlation by figures
windowSize = 50;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

layer= 1;
feature = 39; % 39: middle speed; 13: midbody bend


startpt = 300;% 5300;% 
endpt =  ;%7600; % length(comb_mtx_NaN0{layer}(feature,:)); 

id =  255;

feature_vec = comb_mtx_NaN0{layer}(feature,:);
neuron_vec = comb_mtx_NaN0{layer}(id,:);

rev_fra = (feature_vec<-0);
zero_fra = abs(feature_vec)< 1e-5;
reversals = [];
start_rev = -1;
end_rev = -1;
for jj = 2:length(comb_mtx_NaN0{layer}(feature,:))
   if rev_fra(jj)&& start_rev<0
       start_rev = jj;
   elseif rev_fra(jj)&& (start_rev>0)&&(end_rev<0)
       continue
   elseif ~rev_fra(jj)&& start_rev>0
       reversals = [reversals; start_rev, jj-1];
       start_rev = -1;
       end_rev = -1;
   end
end

reversal_good = [];
for jj = 1:size(reversals,1);
    if reversals(jj,2)-reversals(jj,1)>1.5
        reversal_good = [reversal_good;reversals(jj,:)];
    end
end

ystart = -0.03;
yend = 0.12;
figure, hold on,
basevalue = ystart;
reversals = reversals - startpt;
for ii = 1: size(reversals,1)
    if reversals(ii,1)>0&&(reversals(ii,1)<(endpt-startpt))
    ha = area([reversals(ii,1) reversals(ii,2)], [yend yend],basevalue);
    ha(1).FaceColor = [1 0.8 1];
    ha(1).EdgeColor = [1 0.8 1];
    ha(1).LineStyle = '-';
    ha(1).LineWidth = 0.1;
    ha(1).ShowBaseLine = 'off';
    end
%    plot([reversals(ii,1) reversals(ii,1)], [ystart yend],'y'); 
%    plot([reversals(ii,2) reversals(ii,2)], [ystart yend],'y'); 
end



feature_vec = feature_vec/norm(feature_vec);
feature_vec(zero_fra) = NaN;
%feature_vec = feature_vec-mean(feature_vec);
hold on,plot(feature_vec(startpt:endpt),'r-','LineWidth',2),   % 3050:6700

hold on,
%    plot(comb_mtx_NaN0{1}(43,5000:6000),'m'),

neuron_vec = neuron_vec/norm(neuron_vec);
neuron_vec = neuron_vec + 0.05;
neuron_vec(zero_fra) = NaN;
plot(neuron_vec(startpt:endpt),'k-','LineWidth',2);

plot(xlim, [0 0], 'r-.')


set(gca,'ytick',[])
