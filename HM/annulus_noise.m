function noise = annulus_noise(N, params)
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
%              annulus_index -- binary array which indicates which modes to
%                               generate noise on. sum(annulus_index) = N_A
%                               since noise_params has to get passed
%                               through the nonstiff handle, this is
%                               FLATTENED to a 1D array.
%                               When reshaped to NxN, this is in FFT
%                               ordering
% Output: noise -- noise, same size as q_h
% 
noise_size = params(1); annulus_index = reshape(params(2:end), [N N]); %extract parameters

m_lim = floor((N-1)/2);

% orthant = normrnd(0,1,[m_lim+1,m_lim+1]); %get noise, each mode the same
orthant1 = exp(1i*2*pi*rand(m_lim+1)); %get noise, each mode the same
orthant2 = exp(1i*2*pi*rand(m_lim+1)); %get noise, each mode the same



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

noise = noise_size*(noise.*annulus_index);
end