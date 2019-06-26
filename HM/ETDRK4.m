function [u,t] = ETDRK4(init,T,M,func_noise,func_s,func_nl,params_s,params_nl)
% Implements the Exponential Time Differencing Runge Kutta 4 method 
% (ETDRK4), an explicit one-step, exponential 4th order integrator.
% Designed to solve stiff ODE of the form du/dt = Au + B(u,t)
% 
% Update is of the form (using step size h):
% 
%       a_n = e^(Ah/2)*u_n + Q*B(u_n,t_n)
%       b_n = e^(Ah/2)*u_n + Q*B(a_n,t_n + h/2)
%       c_n = e^(Ah/2)*a_n + Q*[2*B(b_n,t_n + h/2) - B(u_n,t_n)]
%       u_(n+1) = e^(Ah)*u_n + F1*B(u_n,t_n) + 
%                 2*F2*[B(a_n,t_n + h/2)  + B(b_n,t_n + h/2)] + 
%                 F3*B(c_n,t_n + h)
% with 
%       Q = A^(-1) * [e^(Ah/2)-I]
%       F1 = h^(-2) A^(-3) * [-4 - Ah + e^(Ah)*(4-3Ah+(Ah)^2))]
%       F2 = h^(-2) A^(-3) * [2 + Ah + e^(Ah)*(-2+Ah)]
%       F3 = h^(-2) A^(-3) * [-4 - 3Ah - (Ah)^2 - e^(Ah)*(4-Ah)]
% 
% where Q, F1, F2, F3 are precomputed using contour integration. This code
% can ONLY HANDLE DIAGONAL L, since this allows easy contour integration
% (individual contours component wise)
% 
% 
% Args: init -- initial data, at t=0, assumed that there are two initial
%                values to start the multistep method
%        T -- terminal time
%        M -- number of time steps 
%        func_noise -- noise function handle
%        func_s -- RHS of the ODE that is stiff,linear, corresponds to A
%               note, assumes func_s returns THE MATRIX A
%        func_nl -- RHS of the ODE that is nonlinear, corresponds to B
%               note, assumes func_nl returns THE NONLINEAR FUNC EVAL
%        params_s -- parameters to pass into func_s
%        params_nl -- parameters to pass into func_nl
% Output: u -- numerical solution, of size N x M+1
%         t -- array of time steps, of size M+1
% 
h = T/M; %time step

% store each step as a column vector

u = zeros([size(init,1),2]); t = zeros(1,2); %initialize
u(:,1) = init;

A = func_s(u(:,1),params_s); %A only needs to be computed once
A = diag(A); %since system is diagonal, only save this

[Q,F1,F2,F3] = ETDRK4_contour_helper(A,h);
E1 = exp(h*A); E2 = exp(h*A/2); %exponentials we will reuse
disp('contour integral precomps done')

% main solver loop
for i=1:M
    curr_u = u(:,1); %get current solution
    B_init = func_nl(curr_u,func_noise,@nonlinJ,params_nl);
    a = E2.*curr_u + Q.*B_init;
    B_a = func_nl(a,func_noise,@nonlinJ,params_nl);
    b = E2.*curr_u + Q.*B_a;
    B_b = func_nl(b,func_noise,@nonlinJ,params_nl);
    c = E2.*curr_u + Q.*(2*B_b - B_init);
    B_c = func_nl(c,func_noise,@nonlinJ,params_nl);
    
    %update
    u(:,2) = E1.*curr_u + F1.*B_init + 2*F2.*(B_a + B_b) + F3.*B_c;
    
    t(2) = t(1) + h;
    u(:,1) = u(:,2); t(1) = t(2);
end
u = u(:,end); t = t(end);
end



function [Q,F1,F2,F3] = ETDRK4_contour_helper(A,h,n_circ_pts)
% Helper function to use contour integration to get around precomputations
% with removable singularities. Since we assume A is diagonal, each
% component can use a shifted unit circle as the contour, hence the value
% is just the mean around equispaced on the circle
% 
% Args: A -- the diagonal matrix, assume input as a vector of the diagonal
%       h -- step size of the method
%       n_circ_pts -- the number of points to use only the circle for the
%                     contour integration. Default to 32
if nargin < 3
    n_circ_pts = 32; %default
end
equi_pts = exp(1i*2*pi*(0:n_circ_pts-1)/n_circ_pts);

if size(A,1) == 1
    A = A.'; %if A is row vector, change to a column vector
end

% this is an N x n_circ_pts array. Each row is the equispaced points on the
% unit circle, centered at that entry in A. Use row-wise mean
Ah_c = (h*A(:,ones(n_circ_pts,1)))+equi_pts(ones(length(A),1),:);

temp_sq = Ah_c.^2; temp_cube = Ah_c.^3; temp_exp = exp(Ah_c); %temporary variables

%note we use exp for entrywise since diagonal, NOT expm
Q = h*mean( (exp(Ah_c/2)-1)./Ah_c ,2);
F1 = h*mean(-4 - Ah_c + temp_exp.*(4-3*Ah_c+temp_sq) ./ temp_cube,2);
F2 = h*mean(2 + Ah_c + temp_exp.*(-2+Ah_c) ./ temp_cube,2);
F3 = h*mean(-4 - 3*Ah_c - temp_sq + temp_exp.*(4-Ah_c) ./ temp_cube,2);

if norm(imag(A)) == 0 %since if A is real, so are Q,F1,F2,F3
    Q = real(Q); F1 = real(F1); F2 = real(F2); F3 = real(F3);
end
end
