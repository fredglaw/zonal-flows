close all; 
% Script file to compute the energy spectra. We do this both at time slices
% and in statistical equilibrium. For both options we compute the spectra
% both in radial averaged modes and in zonal modes
% 
% To run the simulation we only need to keep running this script. Key flags are:
%       multistep_flag -- 1 is AB2BDF2 (multistep), 0 is ETDRK4
%       real_noise -- 1 is true white noise (SDE solve), 0 is SUF-HM 
%       saver -- 1 means to save q data + vorticity plots, 0 is to not
%       modified -- 1 means mHM, 0 means oHM (with extra Laplacian dissipation)
% 
% Hyperviscosity tuning: N=64   <-->  hype_visc = 7e-23
%                        N=128  <-->  hype_visc = 7e-23 or 5e-25
%                        N=256  <-->  hype_visc = 5e-25
% 
L = 40; %full width of computational box;
sc = L/(2*pi); %scaling factor to go from [-pi,pi] to [-L/2, L/2]
N = 128; %number of nodes in each direction
hype_visc = 7e-23; %hyperviscosity parameter, default 7e-23
gamma = 8; %power on laplacian for hyperviscosity term
kappa = 1; %mean density gradient
alpha = 5; %adiabaticity parameter

x = linspace(-L/2,L/2,N+1); x(end) = []; %delete last entry
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%% Change these from call to call %%%%%%%%%%%%%%%%%%%%
multistep_flag = 0; %flag to see whether to use multistep, AB2BDF2 integrator
% note in reality will want to use white noise, deterministic forcing
% simply mimics the effect of the 2-field model HW
real_noise = 0; %flag to see whether to use white noise, or determinisitic forcing
saver_slice = 0; %flag to see whether or not to save the slice spectra
saver_equil = 0; %flag to see whether or not to save the equilibrium spectra
modified = 1; %flag for oHM or mHM.   0 -> oHM      and      1 -> mHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% INITIAL SETUP / PARAMS FOR THE SIM %%%%%%%%%%%%%%%%%%%
rng(1); %set seed;
init_freq = 1; init_T = 0;

% % small in real space IC
% init_q = (1/500)*(rand(size(X))-(1/2));
% init_q_h = reshape(fft2(init_q), [N*N,1]);

% small in Fourier domain IC
% init_q_h = reshape(fft2(real(ifft2((1/500)*(rand(size(X))-(1/2))))), [N*N,1]);

%BALANCED NORMALIZATION FOR ENERGY
init_q_h = reshape(fft2(real(ifft2((1/500)*(rand(size(X))-(1/2))))), [N*N,1])/N;


zero_mode = init_q_h(1); init_q_h = init_q_h(2:end); %keep zero mode separate
    
if multistep_flag
    init_q_h = [init_q_h, init_q_h];
end

%parameter for the noise size, based on IC, want this independent of current soln to be white in time
if real_noise
	k_f = N/2;
    dk = k_f/8;
    k_vals = -ceil((N-1)/2):floor((N-1)/2); k_sq = k_vals.^2;
    k_full = k_sq + (k_sq');
    annulus_index = (k_full < (k_f + dk)^2) & (k_full > (k_f - dk)^2); %get indices in annulus, FFT ordering
%         eps_param = 1/(2*(k_f^2));
    eps_param = 1/(2*(k_f^2)) * 1e3;
    noise_size = sqrt(2*eps_param*(k_f^2) / (sum(sum(annulus_index))*dt));
    params_noise = [noise_size,reshape(annulus_index,[1, N*N])];
else
    noise_size = (kappa^2)/alpha;
    params_noise = [sc,noise_size];
end

params_ns = [kappa,sc,zero_mode,modified,params_noise];
params_s = [hype_visc,gamma,sc];

%pick noise function
if real_noise
    noise_func = @annulus_noise;
    ab2bdf2 = @AB2BDF2_EM;
    etdrk = @ETDRK4_dcorr;
else
    noise_func = @determ_noise;
    ab2bdf2 = @AB2BDF2;
    etdrk = @ETDRK4;
end

k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %wavenumbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_slices = [1,99,150,250,500]; %time slices
N_slices = 200*T_slices; %N slices
final_T = 1400; %terminal time to end simulation
T = 5; %num step each from T_slices(end) to final_T

unique_k = sqrt(get_unique(N))/sc; unique_k(1) = []; %0th mode does not contribute in energy
n_unique = length(unique_k);
SRAES = zeros(length(T_slices),n_unique); % Sliced Radial. Avg. spectra
SZES = zeros(length(T_slices),floor(N/2)); % Sliced Zonal spectra
ERAES = zeros((final_T - sum(T_slices)) / T, n_unique); % Equil. Radial. Avg. spectra
EZES = zeros((final_T - sum(T_slices)) / T, floor(N/2)); % Equil. Zonal spectra

dt = 1/200;

tic;
%temporal integrator for slices
for i=1:length(T_slices)
    %%% RUN TEMPORAL INTEGRATOR %%%
    if multistep_flag
        %AB2BDF2 solution
        [both_q_h,~] = ab2bdf2(init_q_h,T_slices(i),N_slices(i),noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns);
        both_q_h = [zero_mode zero_mode; both_q_h]; q_h = reshape(both_q_h(:,end),[N,N]); 
%         q = ifft2(q_h); %put zero mode back
        q = ifft2(q_h)*N; %BALANCED NORMALIZATION FOR ENERGY
    else
        %ETDRK2 solution
        [q_h,~] = etdrk(init_q_h,T_slices(i),N_slices(i),noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns);
        q_h = [zero_mode; q_h]; q_h = reshape(q_h,[N,N]); 
%         q = ifft2(q_h); %put zero mode back in
        q = ifft2(q_h)*N; %BALANCED NORMALIZATION FOR ENERGY
    end
    reshape(q_h,[N*N,1]); term_T = init_T + T_slices(i);
    % disp(['Integrating forward by T=',num2str(T),' with M=',num2str(N_time),' steps took t=',num2str(t),' seconds']);
    disp(['Current time we have integrated to is term_T=', num2str(term_T)]);
    disp(['largest imag value is ',num2str(max(max(abs(imag(q)))))]);

    phi_h = - q_h ./ (1 + (k_vals.^2) + ((k_vals').^2));
    if modified
        phi_h(1,2:end) = -q_h(1,2:end) ./ (k_vals(2:end).^2);
        phi_h(1,1) = 0;
    end

    if isnan(max(max(abs(imag(q))))) == 0 %check for blow up
        %update initial and terminal times
        init_T = init_T + T_slices(i); 
        init_q = real(q);

        if multistep_flag
    %         init_q_h = both_q_h(2:end,:);
             init_q_h = [reshape(fft2(real(ifft2(reshape(both_q_h(:,1),[N,N])))),[N*N,1]),...
                 reshape(fft2(real(ifft2(reshape(both_q_h(:,2),[N,N])))),[N*N,1])]/N;
             init_q_h = init_q_h(2:end,:);
        else
            init_q_h = reshape(fft2(real(q)),[N*N,1])/N;
            init_q_h = init_q_h(2:end);
        end
    else
        disp(['got NaNs, try again and integrate a smaller T']);
    end
    %%%%%%
    
    SRAES(i,:) = get_RA_energy(phi_h,sc); %update ZAMFTS
    SZES(i,:) = get_zonal_energy(phi_h,sc); %update TKETS
end


%now in equilibrium, step by T=5 from sum(T_slices) to final_T
N_time = T*200;
for i=1:round((final_T-sum(T_slices))/T)
    %%% RUN TEMPORAL INTEGRATOR %%%
    if multistep_flag
        %AB2BDF2 solution
        [both_q_h,~] = ab2bdf2(init_q_h,T,N_time,noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns);
        both_q_h = [zero_mode zero_mode; both_q_h]; q_h = reshape(both_q_h(:,end),[N,N]); 
%         q = ifft2(q_h); %put zero mode back
        q = ifft2(q_h)*N; %BALANCED NORMALIZATION FOR ENERGY
    else
        %ETDRK2 solution
        [q_h,~] = etdrk(init_q_h,T,N_time,noise_func,@HM_stiff,@HM_nonstiff,params_s,params_ns);
        q_h = [zero_mode; q_h]; q_h = reshape(q_h,[N,N]); 
%         q = ifft2(q_h); %put zero mode back in
        q = ifft2(q_h)*N; %BALANCED NORMALIZATION FOR ENERGY
    end
    reshape(q_h,[N*N,1]); term_T = init_T + T;
    % disp(['Integrating forward by T=',num2str(T),' with M=',num2str(N_time),' steps took t=',num2str(t),' seconds']);
    disp(['Current time we have integrated to is term_T=', num2str(term_T)]);
    disp(['largest imag value is ',num2str(max(max(abs(imag(q)))))]);

    phi_h = - q_h ./ (1 + (k_vals.^2) + ((k_vals').^2));
    if modified
        phi_h(1,2:end) = -q_h(1,2:end) ./ (k_vals(2:end).^2);
        phi_h(1,1) = 0;
    end

    if isnan(max(max(abs(imag(q))))) == 0 %check for blow up
        %update initial and terminal times
        init_T = init_T + T; 
        init_q = real(q);

        if multistep_flag
    %         init_q_h = both_q_h(2:end,:);
             init_q_h = [reshape(fft2(real(ifft2(reshape(both_q_h(:,1),[N,N])))),[N*N,1]),...
                 reshape(fft2(real(ifft2(reshape(both_q_h(:,2),[N,N])))),[N*N,1])]/N;
             init_q_h = init_q_h(2:end,:);
        else
            init_q_h = reshape(fft2(real(q)),[N*N,1])/N;
            init_q_h = init_q_h(2:end);
        end
    else
        disp(['got NaNs, try again and integrate a smaller T']);
    end
    %%%%%%
    
    ERAES(i,:) = get_RA_energy(phi_h,sc); %update ERAES
    EZES(i,:) = get_zonal_energy(phi_h,sc); %update EZES
end
ERAES = mean(ERAES,1); EZES = mean(EZES,1); %mean over equilibrium samples

t=toc;
disp(['Integrating forward to T=',num2str(final_T),' with dt=',num2str(dt),' took t=',num2str(t),' seconds']);


%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%
n_contours = 20;
filled_plot = 1; %flag to do contour filled
colo = parula;
colo = jet;
[T_ts,X_ts] = meshgrid(linspace(T,final_T,round(final_T/T)),x);

figure(1);
for i=1:length(T_slices)
    loglog(unique_k,SRAES(i,:)); hold on;
end
legend('t=1','t=100','t=250','t=500','t=1000')
title('energy spectra, radially averaged')
xlabel('wavenumber');

if saver_slice
    savefig(['SRAES,N',num2str(N),'.fig']);
end

figure(2);
for i=1:length(T_slices)
    loglog((1:(floor(N/2)))/sc,SZES(i,:)); hold on;
end
legend('t=1','t=100','t=250','t=500','t=1000')
title('energy spectra, zonal states')
xlabel('wavenumber');

if saver_slice
    savefig(['SZES,N',num2str(N),'.fig']);
end


figure(3);
loglog(unique_k,ERAES);
title('equilibirum energy spectra, radially averaged')
xlabel('wavenumber');

if saver_equil
    savefig(['ERAES,N',num2str(N),'.fig']);
end

figure(4);
loglog((1:(floor(N/2)))/sc,EZES);
title('equilibirum energy spectra, zonal states')
xlabel('wavenumber');

if saver_equil
    savefig(['EZES,N',num2str(N),'.fig']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%
function KE = get_RA_energy(phi_h,sc)
% Subroutine to get radially averaged energy
N = size(phi_h,1);
k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2));
k_temp = ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %NOT EIGENVALUES, used only for accumarray
k_sq = k_vals.^2 + (k_vals.^2)'; %true k^2
E_vals = k_sq .* (abs(phi_h).^2); %true KE values
k_temp_sq = k_temp.^2 + (k_temp.^2)'; %temp k^2
k_temp_sq = reshape(k_temp_sq,[N*N,1]);  E_vals = reshape(E_vals, [N*N,1]); %reshape for accumarray
KE = [E_vals(1,1); accumarray(k_temp_sq(2:end), E_vals(2:end))];
KE = KE(KE~=0)'; %the correct sorted energies, corresponding to the k values of get_unique(N)
end

function KE = get_zonal_energy(phi_h,sc)
% Subroutine to get zonal states energy
N = size(phi_h,1);
k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2));
k_temp = ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %NOT EIGENVALUES, used only for accumarray
k_sq = k_vals.^2 + (k_vals.^2)'; %true k^2
E_vals = k_sq .* (abs(phi_h).^2); E_vals = E_vals(1,:); %true KE values
k_temp_sq = k_temp.^2; %temp k^2
KE = [E_vals(1); accumarray(k_temp_sq(2:end)', E_vals(2:end))];
KE = KE(KE~=0)'; %the correct sorted energies, corresponding to the k values from 0 to ceil((N-1)/2);
end

function unique_k = get_unique(N)
k_vals = ifftshift(-ceil((N-1)/2):floor((N-1)/2));
k_sq = k_vals.^2 + (k_vals.^2)';
unique_k = unique(k_sq);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


