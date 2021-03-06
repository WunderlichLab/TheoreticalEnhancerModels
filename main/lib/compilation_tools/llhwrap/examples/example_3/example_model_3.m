clear
close all
%clc
%% COMPILATION

[exdir,~,~]=fileparts(which('example_model_3.m'));
% compile the model
llhwrap('model_example_3','example_model_3_syms',exdir)

%% SIMULATION

% time vector
t = [linspace(0,4,21)];
p = [1.1,0.3,1];
k = [];

D.Y = [     1.0171
    1.1761
    1.1680
    1.1359
    1.1778
    1.3423
    1.3079
    1.2784
    1.4976
    1.5903
    1.6585
    1.4688
    1.0999
    1.0128
    0.7198
    0.9814
    0.6755
    0.5091
    0.4471
    0.5249
    0.3288];

D.Sigma_Y = 0.1*ones(size(D.Y));


options.sensi = 1;
options.cvode_maxsteps = 1e6;
options.cvode_rtol = 1e-12;
options.cvode_atol = 1e-12;
% load mex into memory
sol = llh_model_example_3(t,log10(p),k,D,options);

%% Plot

figure
subplot(2,1,1)
errorbar(t,D.Y,D.Sigma_Y)
hold on
% plot(t,sol.y)

xlabel('time t')
ylabel('observable')
title(['log-likelihood: ' num2str(sol.llh) ])

y = (p(2)*t + p(3)).*(t<2) + ( (2*p(2)+p(3)-p(2)/p(1))*exp(-p(1)*(t-2))+p(2)/p(1) ).*(t>=2);


tfine = linspace(0,4,100001);
xfine = (p(2)*tfine + 1).*(tfine<2) + ( (2*p(2)+p(3)-p(2)/p(1))*exp(-p(1)*(tfine-2))+p(2)/p(1) ).*(tfine>=2);

mu = zeros(1,length(tfine));
for it = 1:length(t)
if(t(it)<=2)
mu = mu + ((y(it)-D.Y(it))/(D.Sigma_Y(it)^2))*(tfine<=t(it));
else
mu = mu + ((y(it)-D.Y(it))/(D.Sigma_Y(it)^2))*exp(p(1)*(tfine-t(it))).*(tfine<=t(it)).*(tfine>2) + ((y(it)-D.Y(it))/(D.Sigma_Y(it)^2))*exp(p(1)*(2-t(it))).*(tfine<t(it)).*(tfine<=2);
end
end
plot(tfine,xfine)
legend('data','simulation')
xlim([min(t)-0.5,max(t)+0.5])
subplot(2,1,2)
plot(tfine,mu)
ylabel('adjoint')
xlim([min(t)-0.5,max(t)+0.5])

grad(1,1) = -trapz(tfine,-mu.*xfine.*(tfine>2))*p(1)*log(10);
grad(2,1) = -trapz(tfine,mu)*p(2)*log(10);
grad(3,1) = -mu(1)*p(3)*log(10);

%% FD

eps = 1e-5;
xi = log10(p);
grad_fd_f = NaN(3,1);
grad_fd_b = NaN(3,1);
for ip = 1:3;
    options.sensi = 0;
    xip = xi;
    xip(ip) = xip(ip) + eps;
    solp = llh_model_example_3(t,xip,k,D,options);
    grad_fd_f(ip,1) = (solp.llh-sol.llh)/eps;
    xip = xi;
    xip(ip) = xip(ip) - eps;
    solp = llh_model_example_3(t,xip,k,D,options);
    grad_fd_b(ip,1) = -(solp.llh-sol.llh)/eps;
end

figure
scatter(abs(grad_fd_f),abs(grad))
hold on
scatter(abs(grad_fd_b),abs(grad))
scatter(mean([abs(grad_fd_b),abs(grad_fd_f)],2),abs(grad))
scatter(abs(sol.sllh),abs(grad))
plot([1e1,1e2],[1e1,1e2],'k:')
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('FD_f','FD_b','FD_c','asa')
xlabel('reference absolute gradient value')
ylabel('computed absolute gradient value')

% tic
% for k=1:1000
% options.sensi = 0;
% sol = llh_model_example_3(t,log10(p),k,D,options);
% end
% toc;
% 




