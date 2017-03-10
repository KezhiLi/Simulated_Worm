%%
windowSize = 200;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

layer= 1;
%feature = 65; % 39: middle speed; 13: midbody bend
for ii =1:260;

figure, plot(c_L{layer}(ii,:),'r');   % 3050:6700
hold on,
plot(tp2,'b');

end

y0 = filter(b,a,c_L{layer}(ii,:)); %3050:6700
y = filter(b,a,tp2);%3050:6700

figure, plot(y0,'r'),
hold on,
plot(y)
plot(zeros(1,length(y)))