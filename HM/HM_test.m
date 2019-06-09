close all; clear variables;
% Script file to test running oHM
L = 8*pi; %full width of computational box;
sc = L/(2*pi); %scaling factor to go from [-pi,pi] to [-L/2, L/2]
N = 64; %number of nodes in each direction
hype_visc = 7e-21 ; %hyperviscosity parameter
gamma = 8; %power on laplacian for hyperviscosity term
kappa = 1;
T = 10; %terminal time
M = 10; %number of time steps

x = linspace(-L/2,L/2,N+1); x(end) = []; %delete last entry
[X,Y] = meshgrid(x,x);

init_freq = 1;
%initial condition
init_q = 1*(sin(init_freq*X).*sin(init_freq*Y/sc)); %random initial condition on q
init_q_h = reshape(fft2(init_q), [N*N,1]);
zero_mode = init_q_h(1); init_q_h = init_q_h(2:end); %keep zero mode separate

% build parameters
modified = 1; %flag for oHM or mHM.   0 -> oHM      and      1 -> mHM
params_noise = [sc,1]; %parameters depending on the noise function used
params_ns = [kappa,sc,zero_mode,modified,params_noise];
params_s = [hype_visc,gamma,sc];

%ETDRK2 solution
tic; [q_h,~] = ETDRK2(init_q_h,T,M,@middle_k_noise,@HM_stiff,@HM_nonstiff,params_s,params_ns); t = toc;
q_h = [zero_mode; q_h]; q = ifft2(reshape(q_h,[N,N])); %put zero mode back in


% %RK23 solution
% rhs23s = @(a,u) rhs23s_long(u,@middle_k_noise,@HM_stiff,@HM_nonstiff,params_s,params_ns);
% tic; [q_h,~] = ode23s(rhs23s,linspace(0,T,100),init_q_h); t = toc;
% q_h = [zero_mode; q_h]; q = ifft2(reshape(q_h,[N,N])); %put zero mode back in




disp(['Solving to T=',num2str(T),' with M=',num2str(M),' steps took t=',num2str(t),' seconds']);
disp(['largest imag value is ',num2str(max(max(abs(imag(q)))))]);

%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%
figure(1);
contour(X,Y,real(init_q),25); colorbar;
title('initial q');

figure(2);
contour(X,Y,real(q),25); colorbar;
title(['terminal q at T=',num2str(T)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%
function RHS = rhs23s_long(u,noise_fn,stiff_fn,nonstiff_fn,params_s,params_ns)
% Subroutine to wrap our function handles into the built-in ode23s solver
% 
% Args: u -- current in the ODE solver
%       noise_fn -- noise function handle
%       stiff_fn -- stiff function handle
%       nonstiff_fn -- nonstiff function handle
%       params_s -- parameters for stiff function handle
%       params_ns -- parameters for nonstiff function handle
% Output: RHS -- RHS of the HM equation
% 
A = stiff_fn(u,params_s); %sparse, large stiff matrix
RHS = A*u + nonstiff_fn(u,noise_fn,@nonlinJ,params_ns);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
