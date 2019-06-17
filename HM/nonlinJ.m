function J_h = nonlinJ(phi_h,q_h,sc,aa)
% Routine which returns the product (phi_x)(q_y) - (phi_y)(q_x), in Fourier
% space, using anti-aliasing with 3N/2 rule
% 
% Assume inputs are shaped like N x N, so as 2D arrays
% 
% The ending _h signaling 'hat' meanas representation in Fourier space
% 
% Args: phi_h -- phi in Fourier space
%       q_h -- phi in Fourier space
%       sc -- scaling parameter for the model, takes [-pi,pi] to [-L/2,L/2]
%       aa -- flag for anti-aliasing. Default to true
%                   aa = 1 -> Use anti-aliasing
%                   aa = 0 -> Do not use anti-aliasing
% Output: prod_h -- the product in Fourier space
% 
if nargin < 4
    aa = 1; %by default use antialiasing
end
N = size(phi_h,1); %number of nodes
M = ceil(3*N/2); %finer grid by 3/2 factor
% M = 2*N;

inv_sc = 1/sc; %inverse scaling, jacobian when doing derivatives on [-pi,pi]
% can also be seen as scaling for non-integer wavenumbers/vectors

%the range of wavenumbers, in Fourier ordering
k_vals = -ceil((N-1)/2):floor((N-1)/2);
if aa %when using anti-aliasing, will oversample soln, so 
    if rem(N,2) == 0
        aa_deriv_flag = 1; %flag to note that we are differentiating an unmatched mode
    else
        aa_deriv_flag = 0;
    end
else
    if rem(N,2) == 0
        k_vals(1) = 0;
    end
    aa_deriv_flag = 0;
end
k_vals = ifftshift(k_vals);
% k_vals = inv_sc*ifftshift(k_vals);


% build x,y derivatives, in Fourier ordering
phi_x_h = 1i*k_vals.*phi_h; %column broadcasting for phi_x derivs
phi_y_h = 1i*phi_h.*(k_vals'); %row broadcasting for phi_y derivs
q_x_h = 1i*k_vals.*q_h; %column broadcasting for q_x derivs
q_y_h = 1i*q_h.*(k_vals'); %row broadcasting for q_y derivs

if aa
    J_h = aa_prod(phi_x_h,q_y_h,N,M,aa_deriv_flag,0,1) - aa_prod(phi_y_h,q_x_h,N,M,aa_deriv_flag,1,0); %get J_h
else
    J_h = fft2(ifft2(phi_x_h).*ifft2(q_y_h))-fft2(ifft2(phi_y_h).*ifft2(q_x_h));
end

J_h = (inv_sc^2)*J_h;
end