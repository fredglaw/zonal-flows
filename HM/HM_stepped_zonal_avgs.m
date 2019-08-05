close all; 
% Script file to produce time series of zonally averaged flows for HM.
% Effectively mimics the script HM_stepped_test.m which lets us piecewise
% run the simulation. We begin this from the beginning and hope to see the
% development of the zones through the time series as time goes on.
% 
% To run the simulation we only need to keep running this script. Key flags are:
%       multistep_flag -- 1 is AB2BDF2 (multistep), 0 is ETDRK4
%       real_noise -- 1 is true white noise (SDE solve), 0 is SUF-HM 
%       saver -- 1 means to save q data + vorticity plots, 0 is to not
%       modified -- 1 means mHM, 0 means oHM (with extra Laplacian dissipation)
% 
% Hyperviscosity tuning: N=64   <-->  hype_visc = 7e-23
%                        N=128  <-->  hype_visc = 7e-23
%                        N=256  <-->  hype_visc = 5e-25
% 
L = 40; %full width of computational box;
sc = L/(2*pi); %scaling factor to go from [-pi,pi] to [-L/2, L/2]
N = 256; %number of nodes in each direction
hype_visc = 5e-25; %hyperviscosity parameter, default 7e-23
gamma = 8; %power on laplacian for hyperviscosity term
kappa = 1; %mean density gradient
alpha = 5; %adiabaticity parameter
final_T = 1400;
T = 5; %time increment to simulate with 
N_time = T*200; %number of time steps
dt = T/N_time;

x = linspace(-L/2,L/2,N+1); x(end) = []; %delete last entry
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%% Change these from call to call %%%%%%%%%%%%%%%%%%%%
multistep_flag = 0; %flag to see whether to use multistep, AB2BDF2 integrator
% note in reality will want to use white noise, deterministic forcing
% simply mimics the effect of the 2-field model HW
real_noise = 0; %flag to see whether to use white noise, or determinisitic forcing
saver = 1; %flag to see whether or not to save ZAFTS
modified = 1; %flag for oHM or mHM.   0 -> oHM      and      1 -> mHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% INITIAL SETUP / PARAMS FOR THE SIM %%%%%%%%%%%%%%%%%%%
rng(1); %set seed;
init_freq = 1; init_T = 0;
%initial condition
init_q = (1/500)*(rand(size(X))-(1/2));
init_q_h = reshape(fft2(init_q), [N*N,1]);
zero_mode = init_q_h(1); init_q_h = init_q_h(2:end); %keep zero mode separate
    
if multistep_flag
    init_q_h = [init_q_h, init_q_h];
end

%parameter for the noise size, based on IC, want this independent of current soln to be white in time
if real_noise
	noise_size = 1/sqrt(2*dt);
	noise_size = noise_size*(norm(init_q(:,end))*(L/N))*(1e-8);
else
    noise_size = (kappa^2)/alpha;
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

k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %wavenumbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ZAMFTS = zeros(N,round(final_T/T)); %Zonally Averaged Mean Flow Time Series

init_phi_h = -reshape([zero_mode;init_q_h(:,end)],[N,N]) ./ (1 + (k_vals.^2) + ((k_vals').^2)); % get the vorticity

% ZAMFTS(1,:) = get_zamf(init_phi_h,sc);

tic;
%temporal integrator
for i=1:round(final_T/T)
    %%% RUN TEMPORAL INTEGRATOR %%%
    if multistep_flag
        %AB2BDF2 solution
        [both_q_h,~] = ab2bdf2(init_q_h,T,N_time,noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns);
        both_q_h = [zero_mode zero_mode; both_q_h]; q_h = reshape(both_q_h(:,end),[N,N]); q = ifft2(q_h); %put zero mode back
    else
        %ETDRK2 solution
        [q_h,~] = etdrk(init_q_h,T,N_time,noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns);
        q_h = [zero_mode; q_h]; q_h = reshape(q_h,[N,N]); q = ifft2(q_h); %put zero mode back in
    end
    reshape(q_h,[N*N,1]); term_T = init_T + T;
    % disp(['Integrating forward by T=',num2str(T),' with M=',num2str(N_time),' steps took t=',num2str(t),' seconds']);
    disp(['Current time we have integrated to is term_T=', num2str(term_T)]);
    disp(['largest imag value is ',num2str(max(max(abs(imag(q)))))]);

    phi_h = - q_h ./ (1 + (k_vals.^2) + ((k_vals').^2));

    if isnan(max(max(abs(imag(q))))) == 0 %check for blow up
        %update initial and terminal times
        init_T = init_T + T; 
        init_q = real(q);

        if multistep_flag
    %         init_q_h = both_q_h(2:end,:);
             init_q_h = [reshape(fft2(real(ifft2(reshape(both_q_h(:,1),[N,N])))),[N*N,1]),...
                 reshape(fft2(real(ifft2(reshape(both_q_h(:,2),[N,N])))),[N*N,1])];
             init_q_h = init_q_h(2:end,:);
        else
            init_q_h = reshape(fft2(real(q)),[N*N,1]);
            init_q_h = init_q_h(2:end);
        end
    else
        disp(['got NaNs, try again and integrate a smaller T']);
    end
    %%%%%%
    
    ZAMFTS(:,i) = get_zamf(phi_h,sc);%update ZAMFTS

end
t=toc;
disp(['Integrating forward to T=',num2str(final_T),' with dt=',num2str(dt),' took t=',num2str(t),' seconds']);


%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%
n_contours = 20;
filled_plot = 1; %flag to do contour filled
colo = parula;
colo = jet;
[T_ts,X_ts] = meshgrid(linspace(T,final_T,round(final_T/T)),x);

figure(1);
contourf(T_ts,X_ts,real(ZAMFTS),n_contours); colorbar;
title(['ZAMFTS at T=',num2str(init_T)]);


if saver
    savefig(['ZAFTS,T',num2str(term_T),',N',num2str(N),'.fig']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%
function zamf = get_zamf(phi_h,sc)
% Subroutine to zonally average the current flow given the electrostatic
% potential in Fourier space.
N = size(phi_h,1);
k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2));
zamf = (ifft(1i*k_vals.*phi_h(1,:))/N)';
% zamf = (ifft(phi_h(1,:))/N)';

% top row of phi_h encodes the zonally averaged ES potential. We spectrally
% differentiate, and divide by N to account for ifft instead of ifft2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


