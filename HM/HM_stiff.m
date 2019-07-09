function RHS_stiff_h = HM_stiff(q_h,params)
% Routine to return the nonlinear piece of RHS of the oHM with a
% noise function handle, and a nonlinear interaction function handle. 
% Manually enforce the consistency condition that q is the potential vorticity
% 
% Returns the matrix for the system, sparse
% 
% Assume q_h comes in the shape of N^2 x 1
% 
% Args: q_h -- q in Fourier space
%       params -- vector of parameters. Contains:
%                  hype_visc -- hyperviscosity parameter
%                  gamma -- power of Laplacian parameter
%                  sc -- scaling parameter for the model, takes [-pi,pi] to [-L/2,L/2]
% Output: RHS_stiff_h -- RHS of oHM PDE, stiff piece, Fourier space
% 
hype_visc = params(1); gamma = params(2); sc = params(3); %extract parameters

N = round(sqrt(length(q_h) + 1)); %number of nodes
k_vals = (1/sc)*ifftshift(-ceil((N-1)/2):floor((N-1)/2)); %wavenumbers
k_vals_sq = k_vals.^(2);

coeffs = hype_visc*((-1)^gamma)*((k_vals_sq + (k_vals_sq')).^gamma); %get stiff coeffs

%reshape in N^2 x 1 size, then make the diagonal of an N^2 x N^2 matrix
temp_coeffs = reshape(coeffs,[N*N,1]);
RHS_stiff_h = spdiags(temp_coeffs(2:end),0,N*N - 1,N*N - 1); %do not include first term (0,0) mode
end

