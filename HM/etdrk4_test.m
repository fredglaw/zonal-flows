% script to test the etdrk4 accuracy
% do it on the standard du/dt = Au, with no nonlinear term
A = 1;
T = 5; %terminal time

blank_noise = @(u,params) zeros(size(u));
nonlin = @(u,noise,excess,params) zeros(size(u));
lin = @(u,params) A*ones(length(u));

init = 1;



dts = 2.^(-(0:3));
Ms = T./dts;
errs = zeros(size(dts));
ms_errs = zeros(size(errs));
truth = exp(A*T);
for i=1:length(dts)
    errs(i) = abs(ETDRK4(init,T,Ms(i),blank_noise,lin,nonlin,[],[])-truth)./abs(truth);
    temp = AB2BDF2([exp(-dts(i)),1],T,Ms(i),blank_noise,lin,nonlin,[],[]);
    ms_errs(i) = abs(temp(:,end)-truth)./abs(truth);
end