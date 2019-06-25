function RHS_nonstiff_h = HM_nonstiff(q_h,noise_fn,nonlin_fn,params)
% Routine to return the nonstiff piece of RHS of the oHM with a
% noise function handle, and a nonstiff interaction function handle. 
% Manually enforce the consistency condition that q is the potential vorticity
% 
% Since this should be fed into some temporal integrator, while otherwise
% we have been storing q, phi as 2D arrays, we flatten to an (N^2) x 1
% array in the output
% 
% Hence also assume q_h comes in the shape of N^2 x 1
% 
% Args: q_h -- q in Fourier space
%       noise_fn -- some input noise function handle, in Fourier space
%       nonlin_fn -- function handle for nonlinear interaction of q and 
%                    phi, in Fourier space
%       params -- vector of parameters. Contains:
%                  kappa -- parameter for y-deriv of phi
%                  sc -- scaling parameter for the model, takes [-pi,pi] to [-L/2,L/2]
%                  zero_mode -- the zeroth mode of th IC
%                  noise_params -- params for the noise function
%                  modified -- flag for which elliptic PDE links q and phi.
%                                      modified=1 -> modified HM
%                                      modified=0 -> original HM
% Output: RHS_nonstiff_h -- RHS of oHM PDE, nonstiff piece, Fourier space
%                         Note this is output as (N^2) x 1 vector
% 
kappa = params(1); sc = params(2); zero_mode = params(3); modified = params(4);
noise_params = params(5:end); %extract params

q_h = [zero_mode;q_h]; %add zero mode back in for the computation
N = round(sqrt(length(q_h))); %number of nodes
q_h = reshape(q_h,[N,N]); %reshape
k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %wavenumbers
full_ks = ((k_vals.^2) + ((k_vals.^2)'));
phi_h = -q_h ./ (1 + full_ks); % consistency, q is potential vorticity

% if mHM instead of oHM
if modified
    phi_h(1,2:end) = -q_h(1,2:end) ./ (k_vals(2:end).^2);
    phi_h(1,1) = 0;
end

J_h = nonlin_fn(phi_h,q_h,sc);
if rem(N,2) == 0
    k_vals = fftshift(k_vals); k_vals(1) = 0; k_vals = ifftshift(k_vals);
end
phi_y_h = (1i*phi_h.*(k_vals')); %left broadcasting for y-derivs


% the right hand side, before noise
RHS = -J_h + kappa*phi_y_h;

%%% ONLY IF FLUX BALANCED
RHS = RHS - (5e-4)*full_ks.*q_h;

%%% ONLY IF DETERMINISTIC NOISE
noise = noise_fn(q_h,noise_params); %generate the noise
RHS = RHS + noise;


% disp(['norm imag value of RHS is ',num2str(norm(imag(ifft2(RHS))))]);

% % go back to real space to enforce real
% RHS = fft2(real(ifft2(RHS)));

reshaped_RHS = reshape(RHS,[N*N, 1]);
RHS_nonstiff_h = reshaped_RHS(2:end);
% RHS_nonstiff_h = reshape(RHS(2:end), [N*N - 1,1]);
end

