%%
windowSize = 200;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

layer= 1;
feature = 39; % 39: middle speed; 13: midbody bend
for ii = 56+1:5:56+260;
close

scrsz = get(groot,'ScreenSize');
figure('Position',[10 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5])

fold = sqrt(norm(comb_mtx_NaN0{layer}(feature,:)));
subplot(6,1,1)
plot(comb_mtx_NaN0{layer}(ii,1:5000)*fold,'r');   % 3050:3200
hold on, plot(comb_mtx_NaN0{layer}(feature,1:5000),'b');   % 3050:3200
subplot(6,1,2)
plot(comb_mtx_NaN0{layer}(ii+1,1:5000)*fold,'r');   % 3050:3200
hold on, plot(comb_mtx_NaN0{layer}(feature,1:5000),'b');   % 3050:3200
subplot(6,1,3)
plot(comb_mtx_NaN0{layer}(ii+2,1:5000)*fold,'r');   % 3050:3200
hold on, plot(comb_mtx_NaN0{layer}(feature,1:5000),'b');   % 3050:3200
subplot(6,1,4)
plot(comb_mtx_NaN0{layer}(ii+3,1:5000)*fold,'r');   % 3050:3200
hold on, plot(comb_mtx_NaN0{layer}(feature,1:5000),'b');   % 3050:3200
subplot(6,1,5)
plot(comb_mtx_NaN0{layer}(ii+4,1:5000)*fold,'r');   % 3050:3200
hold on, plot(comb_mtx_NaN0{layer}(feature,1:5000),'b');   % 3050:3200

subplot(6,1,6)
plot(comb_mtx_NaN0{layer}(feature,1:5000),'b');   % 3050:3200
hold off



% plot(res_L{layer}(ii+1,3050:3200),'m');
% plot(res_L{layer}(ii+2,3050:3200),'c');

%plot(tp1,'b');

% figure, plot(res_L{layer}(ii+1,3050:3200),'r');   % 3050:6700
% hold on,
% plot(tp1(3050:3200),'b');
% 
% figure, plot(res_L{layer}(ii+2,3050:3200),'r');   % 3050:6700
% hold on,
% plot(tp1(3050:3200),'b');

end

y0 = filter(b,a,c_L{layer}(ii,:)); %3050:6700
y = filter(b,a,tp2);%3050:6700

figure, plot(y0,'r'),
hold on,
plot(y)
plot(zeros(1,length(y)))