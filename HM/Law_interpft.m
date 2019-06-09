function output = Law_interpft(f,N,deriv,int)
% CODE BY FREDERICK LAW
% 
% INPUT: 'f' -- function to take the FFT (discrete values)
%        'N' -- size of finer grid for evaluation of interpolant,
%        'deriv' -- optional argument determining whether to compute
%                   an approximation to derivative.
%        'int' -- optional argument determining whether to compute
%                   an approximation to integral
% 
n = length(f);
% FFT to get Fourier coefficients, shifted to natural ordering
f_hat = fftshift(fft(f));

% check to see if we are computing the interpolant or the
% DERIVATIVE of the interpolant.
% for derivative, map Fourier coefficients by: f_k  -->  ikf_k
if nargin == 3
    freqs= -ceil((n-1)/2):floor((n-1)/2);
    f_hat = i*freqs.*f_hat;
end

% check to see if we are computing the interpolant or the
% INTEGRAL of the interpolant.
% for integral map Fourier coefficients by: f_k  -->  -if_k / k
% which hold for k nonzero. set f_0 = 0;
if nargin == 4
    f0 = ifftshift(f_hat); f0 = f0(1); %shift, save 0th mode
    freqs= -ceil((n-1)/2):floor((n-1)/2);
    f_hat = -i*f_hat ./ freqs; %rescale modes
    f_hat = ifftshift(f_hat); f_hat(1) = 0;%shift back, set 0th to 0
    f_hat = fftshift(f_hat); %shift back again
end



% interpolate by performing an IFFT, set high frequencies to 0;
% this is where the shift comes in, we can manually embedd f_hat;
% largest negative Fourier coeff corresponds to front_index;
% largest positive Fourier coeff corresponds to back_index;
new_hat = zeros(1, N);


% definition of these placements depends on parity of n
par_n = rem(n,2);
if par_n == 0
    f_hat(1) = f_hat(1)/2; %splitting unmatched mode
    f_hat(end+1) = conj(f_hat(1)); %manually matching the unmatched mode
    front_index = ceil((N-1-n)/2) + 1;
    back_index = front_index+n;
else
    front_index = ceil((N-n)/2) + 1; 
    back_index = front_index+n-1;
end

% zero padding by inserting f_hat into a longer vector
new_hat(front_index:back_index) = f_hat;


% MATLAB implements the scaling in IFFT, so we have to multiply by N
% and then divide by 1/n to get the correct scaling factor
output = (N/n)*ifft(ifftshift(new_hat));

% In the case that we compute the integral, need to manually add in
% the other terms
if nargin == 4
    val = linspace(-pi, pi, N+1); val(end) = [];
    % adding in extra sum + f_0 (x)
    output = output + (1/n)*(-real(sum(f_hat)) + f0*val);
end

end