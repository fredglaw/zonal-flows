close all; 
% Script file to test running oHM
L = 40; %full width of computational box;
sc = L/(2*pi); %scaling factor to go from [-pi,pi] to [-L/2, L/2]
N = 64; %number of nodes in each direction
hype_visc = 7e-21; %hyperviscosity parameter, default 7e-21
gamma = 2*8; %power on laplacian for hyperviscosity term
kappa = 1; %mean density gradient
alpha = 5; %adiabaticity parameter
T = 100; %terminal time
N_time = T*200; %number of time steps
dt = T/N_time;

x = linspace(-L/2,L/2,N+1); x(end) = []; %delete last entry
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%% Change these from call to call %%%%%%%%%%%%%%%%%%%%
is_first_time = 1;
multistep_flag = 0; %flag to see whether to use multistep, AB2BDF2 integrator
real_noise = 0; %flag to see whether to use white noise, or determinisitic forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if is_first_time
    rng(1); %set seed;
    init_freq = 1;
    %initial condition
%     init_q = (1)*(sin(init_freq*X).*sin(init_freq*Y/sc)); %random initial condition on q
%     init_q = (1)*(sin(init_freq*X).*sin(init_freq*Y/sc)) + (1/50)*(rand(size(X))-(1/2)); %random initial condition on q
%     init_q = 1*ones(size(X));
    init_q = (1/500)*(rand(size(X))-(1/2));
    init_q_h = reshape(fft2(init_q), [N*N,1]);
    zero_mode = init_q_h(1); init_q_h = init_q_h(2:end); %keep zero mode separate
    
    init_T = 0;
    
    %in the case of multistep, use the a perturbation on the initial state
    %as the prior position.
    if multistep_flag
        init_q_h = [init_q_h, init_q_h];
    end
end

term_T = init_T + T;


% build parameters
modified = 1; %flag for oHM or mHM.   0 -> oHM      and      1 -> mHM

%parameter for the noise size, based on IC, want this independent of current soln to be white in time
if is_first_time
    if real_noise
        noise_size = 1/sqrt(2*dt);
        noise_size = noise_size*(norm(init_q(:,end))*(L/N))*(1e-8);
    else
        noise_size = (kappa^2)/alpha;
    end
end


params_noise = [sc,noise_size]; %parameters depending on the noise function used
params_ns = [kappa,sc,zero_mode,modified,params_noise];
params_s = [hype_visc,gamma,sc];

%pick noise function
if real_noise
    noise_func = @middle_k_noise;
    ab2bdf2 = @AB2BDF2_EM;
    etdrk = @ETDRK4_EM;
else
    noise_func = @determ_noise;
    ab2bdf2 = @AB2BDF2;
    etdrk = @ETDRK4;
end

%temporal integrator
if multistep_flag
    %AB2BDF2 solution
    tic; [both_q_h,~] = ab2bdf2(init_q_h,T,N_time,noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns); t = toc;
    both_q_h = [zero_mode zero_mode; both_q_h]; q_h = reshape(both_q_h(:,end),[N,N]); q = ifft2(q_h); %put zero mode back
else
    %ETDRK2 solution
    tic; [q_h,~] = etdrk(init_q_h,T,N_time,noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns); t = toc;
    q_h = [zero_mode; q_h]; q_h = reshape(q_h,[N,N]); q = ifft2(q_h); %put zero mode back in
end
reshape(q_h,[N*N,1]);
disp(['Integrating forward by T=',num2str(T),' with M=',num2str(N_time),' steps took t=',num2str(t),' seconds']);
disp(['Current time we have integrated to is term_T=', num2str(term_T)]);
disp(['largest imag value is ',num2str(max(max(abs(imag(q)))))]);


k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %wavenumbers

init_phi_h = -reshape([zero_mode;init_q_h(:,end)],[N,N]) ./ (1 + (k_vals.^2) + ((k_vals').^2)); % get the vorticity
init_zeta = ifft2(-((k_vals.^2) + ((k_vals').^2)).*init_phi_h);

phi_h = - q_h ./ (1 + (k_vals.^2) + ((k_vals').^2));
zeta = ifft2(-((k_vals.^2) + ((k_vals').^2)).*phi_h);



%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%
n_contours = 20;
filled_plot = 1; %flag to do contour filled
colo = parula;
colo = jet;


figure(1); colormap(colo);
if filled_plot
    contourf(X,Y,real(init_q),n_contours);  colorbar;
else
    contour(X,Y,real(init_q),n_contours); colorbar;
end
title(['q at T=',num2str(init_T)]);


figure(2); colormap(colo);
if filled_plot
    contourf(X,Y,real(q),n_contours); colorbar;
else
    contour(X,Y,real(q),n_contours); colorbar;
end
title(['q at T=',num2str(term_T)]);
savefig(['full,etdrk4,T',num2str(term_T),'.fig']);
save(['full,etdrk4,T',num2str(term_T),'.mat'],'q');


figure(3); colormap(colo);
if filled_plot
    contourf(X,Y,real(init_zeta),n_contours); colorbar;
else
    contour(X,Y,real(init_zeta),n_contours); colorbar;
end
title(['vorticity at T=',num2str(init_T)]);


figure(4); colormap(colo);
if filled_plot
    contourf(X,Y,real(zeta),n_contours); colorbar;
else
    contour(X,Y,real(zeta),n_contours); colorbar;
end
title(['vorticity at T=',num2str(term_T)]);




if isnan(max(max(abs(imag(q))))) == 0
    %update initial and terminal times
    init_T = init_T + T; 
    init_q = real(q);
    
    if multistep_flag
%         init_q_h = both_q_h(2:end,:);
         init_q_h = [reshape(fft2(real(ifft2(reshape(both_q_h(:,1),[N,N])))),[N*N,1]),...
             reshape(fft2(real(ifft2(reshape(both_q_h(:,2),[N,N])))),[N*N,1])];
         init_q_h = init_q_h(2:end,:);
    else
        init_q_h = reshape(fft2(q),[N*N,1]);
        init_q_h = init_q_h(2:end);
    end
else
    disp(['got NaNs, try again and integrate a smaller T']);
end
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
%       params_ns -reshape(q_h,[N*N,1]); init_q_h = init_q_h(2:end); - parameters for nonstiff function handle
% Output: RHS -- RHS of the HM equation
% 
A = stiff_fn(u,params_s); %sparse, large stiff matrix
RHS = A*u + nonstiff_fn(u,noise_fn,@nonlinJ,params_ns);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
