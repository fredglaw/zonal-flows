function noise = middle_k_noise(N, params)
% Generates noise for dynamics on the variable q_h in the Fourier domain.
% We assume q_h is 2D array, the noise is i.i.d. Gaussian, then each entry
% is weighted by the pdf of a Gaussian in Fourier space, centered at the
% median wavenumber(wavevector) of q_h.
% 
% Args: N -- the number of nodes in one direction
%       params -- vector of parameters. Contains:
%              sc -- scaling parameter for the model, takes [-pi,pi] to [-L/2,L/2]
%              noise_size -- size of noise, based on IC, indep of current
%                            q_h so as to be white in time
% Output: noise -- noise, same size as q_h
% 
sc = params(1); noise_size = params(2); %extract parameters

m_lim = floor((N-1)/2);
k_vals = (1/sc)*(-ceil((N-1)/2):floor((N-1)/2));
k_vals = ifftshift(k_vals);
m = (1/sc)*floor((N-1)/4); %median of noise, in positive part


% orthant = normrnd(0,1,[m_lim+1,m_lim+1]); %get noise, each mode the same
orthant1 = (1/sqrt(2))*(normrnd(0,1,[m_lim+1,m_lim+1]) + 1i*(normrnd(0,1,[m_lim+1,m_lim+1]))); %get noise, each mode the same
orthant2 = (1/sqrt(2))*(normrnd(0,1,[m_lim+1,m_lim+1]) + 1i*(normrnd(0,1,[m_lim+1,m_lim+1]))); %get noise, each mode the same



%scale modes by pdf of 2D gaussian
one_d_pdf = normpdf(k_vals(1:m_lim+1),m,m/6);
% orthant1 = orthant1 .* (one_d_pdf.* (one_d_pdf')); 
% orthant2 = orthant2 .* (one_d_pdf.* (one_d_pdf')); 


%use symmetry to reflect the noise into each other quadrant
temp1 = [flip(conj(orthant1(2:end,1))),flip(orthant2(2:end,2:end),1)];
% temp1 = flip(orthant2(2:end,:),1); temp1(:,1) = flip(conj(orthant1(2:end,1)));
temp2 = [flip(conj(orthant1(1,2:end)));flip(conj(orthant2(2:end,2:end)),2)];
% noise = [flip(flip(conj(orthant1(2:end,2:end)),1),2), temp1;...
%                                      flip(conj(orthant(:,2:end)),2), orthant1;];
noise = [flip(flip(conj(orthant1(2:end,2:end)),1),2), temp1;...
                                     temp2, orthant1;];
if rem(N,2) == 0
    noise = padarray(noise, [1 1], 0, 'pre');
end
noise = ifftshift(noise); %shift back to FFT ordering
% is_even = (rem(N,2) == 0);
% placeholder_noise = ifftshift([flip(flip(orthant(2:end,2:end),1),2), flip(orthant(2:end,:),1);...
%                                      flip(orthant(:,2:end),2), orthant;]);
% noise(1+is_even:end,1+is_even:end) = placeholder_noise;

% noise = norm(q_h).*ifftshift(noise); %shift back, and scale based on data
noise = noise_size*noise;
end