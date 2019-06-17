function prod_h = aa_prod(term1_h,term2_h,N,M,aa_deriv_flag,term1_axis,term2_axis)
% Subroutine to compute an antialiased product of term1, term2 using M
% terms, assumed unshifted
% 
% Args: term1_h -- first term, in Fourier space
%       term2_h -- second term, in Fourier space
%       N -- original number of nodes
%       M -- number of nodes to oversample with
%       aa_deriv_flag -- optional flag noting that term1 and term2 are
%                        derivatives, and we have an even number of modes, 
%                        so handle the unmatched in oversample
%       term1_axis -- axis along which differentiation happens for term1
%       term2_axis -- axis along which differentiation happens for term2
% Output: prod_h -- the anti-aliased product, in Fourier space
% 
if nargin < 5
    aa_deriv_flag = 0;
    term1_axis = 0; term2_axis = 0;
end
m_lim = floor((N-1)/2); %index limits of the matched modes

% oversample first in x, then in y

% term1 = ifft2(term1_h); term2 = ifft2(term2_h); %get values in real space
% term1_os = interpft(interpft(real(term1),M,1),M,2);
% term2_os = interpft(interpft(real(term2),M,1),M,2);

% % manual padding, faster
term1_os = manual2Doversample(term1_h,N,M,aa_deriv_flag,term1_axis);
term2_os = manual2Doversample(term2_h,N,M,aa_deriv_flag,term2_axis);

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