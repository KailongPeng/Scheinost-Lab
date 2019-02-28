function [modelfun,beta] = kailong_nonlinear(x,y)

modelfun = @(b,x)(b(1)+b(2)*x.^2+b(3)*x.^3+b(4)*x.^4);
rng('default') % for reproducibility
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';

beta0 = [1;1;1;1];
[beta,~,~,~,~,~] = nlinfit(x,y,modelfun,beta0,opts);

% figure;hold on;
% plot(x,y,'.')
% x_temp = [min(x):(max(x)-min(x))/1000:max(x)];
% y_predicted = modelfun(beta,x_temp);
% plot(x_temp,y_predicted);