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
if rem(N,2) == 0
    k_vals(1) = 0; %zero unmatched mode since we will use for derivatives.
end
k_vals = ifftshift(k_vals);
% k_vals = inv_sc*ifftshift(k_vals);


% build x,y derivatives, in Fourier ordering
phi_x_h = 1i*k_vals.*phi_h; %column broadcasting for phi_x derivs
phi_y_h = 1i*phi_h.*(k_vals'); %row broadcasting for phi_y derivs
q_x_h = 1i*k_vals.*q_h; %column broadcasting for q_x derivs
q_y_h = 1i*q_h.*(k_vals'); %row broadcasting for q_y derivs

if aa
    J_h = aa_prod(phi_x_h,q_y_h,N,M) - aa_prod(phi_y_h,q_x_h,N,M); %get J_h
else
    J_h = fft2(ifft2(phi_x_h).*ifft2(q_y_h))-fft2(ifft2(phi_y_h).*ifft2(q_x_h));
end

J_h = (inv_sc^2)*J_h;
end


function prod_h = aa_prod(term1_h,term2_h,N,M)
% Subroutine to compute an antialiased product of term1, term2 using M
% terms, assumed unshifted
% 
% Args: term1_h -- first term, in Fourier space
%       term2_h -- second term, in Fourier space
%       N -- original number of nodes
%       M -- number of nodes to oversample with
% Output: prod_h -- the anti-aliased product, in Fourier space
% 
m_lim = floor((N-1)/2); %index limits of the matched modes

% oversample first in x, then in y

% term1 = ifft2(term1_h); term2 = ifft2(term2_h); %get values in real space
% term1_os = interpft(interpft(real(term1),M,1),M,2);
% term2_os = interpft(interpft(real(term2),M,1),M,2);

% % manual padding, faster
term1_os = manual2Doversample(term1_h,N,M);
term2_os = manual2Doversample(term2_h,N,M);

%since our data is real, we enforce real
prod_os = real(term1_os) .* real(term2_os); %product, using oversampled, in real space
prod_os_h = ((N/M)^2)*fft2(prod_os); %product in Fourier space, oversampled

% extract the correct indices
% prod_h = zeros(size(term1_h));
% prod_h(1:m_lim+1,1:m_lim+1) = prod_os_h(1:m_lim+1, 1:m_lim+1);
% prod_h(m_lim+2:end,m_lim+2:end) = prod_os_h(end-(N-(m_lim+2)):end,end-(N-(m_lim+2)):end);

% get_corner = ceil((M-N)/2);
get_corner = floor((M-N-1)/2) + 1;
if rem(M,2) == 0
    get_corner = get_corner + 1; 
end

prod_os_h = fftshift(prod_os_h);

%now truncate off the higher order modes. in the case of even we need to
%manually condense info BACK to the unmatched modes
if rem(N,2) == 0
    full_trunc = prod_os_h(get_corner:get_corner+N, get_corner:get_corner+N); %odd of size N+1, all modes full represented
    correct_trunc = full_trunc(1:end-1, 1:end-1); %will need to recondense info for the unmatched mode, manually

    correct_trunc(1,1) = full_trunc(1,1) + full_trunc(1,end) + full_trunc(end,1) + full_trunc(end,end); %full unmatched
    correct_trunc(1,2:end) = full_trunc(1,2:end-1) + full_trunc(end,2:end-1); %unmatched in y, matched in x
    correct_trunc(2:end,1) = full_trunc(2:end-1,1) + full_trunc(2:end-1, end); %unmatched in x, matched in y
else %odd case we do not have unmatched modes
    correct_trunc = prod_os_h(get_corner:get_corner+N-1, get_corner:get_corner+N-1);
end

prod_h = ifftshift(correct_trunc);
% prod_h = ifft2(real(fft2(prod_h)));
% disp(norm(imag(ifft2(prod_h))));
end