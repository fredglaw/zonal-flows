function noise = determ_noise(q_h, params)
% Generates "noise" for dynamics on the variable q_h in the Fourier domain.
% We assume q_h is 2D array, the noise is i.i.d. Gaussian, then each entry
% is weighted by the pdf of a Gaussian in Fourier space, centered at the
% median wavenumber(wavevector) of q_h.
% 
% NOTE: if we choose this form of noise, this is to be discretized FULLY as
% part of the RHS, i.e. we do not handle the noise as an EM term for an SDE
% 
% Args: N -- the number of nodes in one direction
%       params -- vector of parameters. Contains:
%              sc -- scaling parameter for the model, takes [-pi,pi] to [-L/2,L/2]
%              noise_size -- size of noise, based on IC, indep of current
%                            q_h so as to be white in time
% Output: noise -- noise, same size as q_h
% 
sc = params(1); noise_size = params(2); %extract parameters

N = size(q_h,1);
k_vals = (1/sc)*(-ceil((N-1)/2):floor((N-1)/2));
k_vals = ifftshift(k_vals);
full_ks = (k_vals.^2) + (k_vals.^2)';

noise = ((full_ks.*((k_vals.^2)')).*q_h)./((full_ks + 1).^3);

noise = noise_size*noise;
end