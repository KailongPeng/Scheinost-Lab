function y = playwithatanh(x)

y =.5*log((1+x)/(1-x));

close all
x = [-1:0.01:1];
y = atanh(x);
plot(x,y)
figure;
y = atanh(y);
plot(x,y)
figure;
y = atanh(y);
plot(x,y)


% x = [-1:0.01:1];
% y = tanh(x);
% plot(x,y)


zrho = atanh(rho);