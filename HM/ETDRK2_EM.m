function [u,t] = ETDRK2_EM(init,T,M,func_noise,func_s,func_nl,params_s,params_nl)
% Implements the Exponential Time Differencing Runge Kutta 2 method 
% (ETDRK2), an explicit one-step, exponential 2nd order integrator, with
% Euler-Maruyama(EM) for the noise term
% 
%  u_(n+1,*) = exp(A*dt) u_n + F1*B(u_n)              <-- predictor
%  u_(n+1) = u_(n+1,*) + F2*[B(u_(n+1),*)- B(u_n)]    <-- corrector
%
% where     F1 = A^(-1) * (exp(Adt) - I)
%           F2 = A^(-2) * (exp(Adt) - I - Adt) / dt
% 
% Assumes ODE is time-homogenous of the form u' = Au + B(u) where
% --A is stiff but linear
% --B is nonlinear
% 
% INPUT: init -- initial data, at t=0, assumed that there are two initial
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
% OUTPUT: u -- numerical solution, of size N x M+1
%         t -- array of time steps, of size M+1
% 
dt = T/M; %time step

% store each step as a column vector
u = zeros([size(init,1),2]); t = zeros(1,2); %initialize
u(:,1) = init;

A = func_s(u(:,1),params_s); %A only needs to be computed once
tol = 1e-10; %tolerance for computing these

%Note to compute F1, F2 in a stable manner, depends on the size. From
%numerical experiments, seems like threshold is if norm(A*dt) > 13.5*pi. If
%so, use direct. Otherwise, use Taylor poly.

% Note this threshold is not sharp, generally works but sometimes might run
% into stability. So solve once, check if there are NaNs, and if so resolve
% with other method.

% flag is variable indicating which method to use in expm_1,2sing
% default is flag = 2, which means choose preferred method
% if that fails, the new flag says which method failed, resolve with other
init_flag = 2;
[test_u,test_t,flag] = solve_once(u,t,M,A,dt,tol,func_noise,func_nl,params_nl,init_flag);

if max(max(isnan(test_u))) == 1 %indicates we have NaN's, redo computation, with flipped flag
    disp('ya goofed, run 2');
    [u,t,~] = solve_once(u,t,M,A,dt,tol,func_noise,func_nl,params_nl,~flag); %resolve, force other method
else
    u = test_u; t = test_t; %accept
end
u = u(:,end); t = t(end);
end


function [u,t,flag] = solve_once(u,t,M,A,dt,tol,func_noise,func_nl,params_nl,flag)
%one solve of the ETDRK2 method
eAdt = expm(A*dt); %compute the full matrix exponentional only once
[F1,~] = expm_1sing(A,dt,tol,flag); %computed w input flag
[F2,flag] = expm_2sing(A,dt,tol,flag); %update flag
disp('done with computing F1, F2');
N = ceil(sqrt(size(u,1))); noise_params = params_nl(5:end);
blank_noise = @(M,y) zeros(M); %empty noise function, need to feed into nonlin

if isdiag(A) %exp(Adt), F1, F2, all diag so  mat-vecs mults become vec-vec mults.
disp('is diag');
eAdt = diag(eAdt);
F1 = diag(F1); F2 = diag(F2);
	for i=1:M
        B1 = func_nl(u(:,1),blank_noise,@nonlinJ, params_nl); %first nonlinear func eval
        predictor = eAdt.*u(:,1) + F1.*B1; %predictor

        %second nonlinear func eval and corrector
        u(:,2) = predictor + F2.*(func_nl(predictor,blank_noise,@nonlinJ, params_nl)-B1);
        
        %EM step
        curr_noise = reshape(func_noise(N,noise_params),[N*N,1]);
        u(:,2) = u(:,2) + dt*curr_noise(2:end);
        
        t(2) = t(1) + dt;
        u(:,1) = u(:,2); t(1) = t(2);
	end
else %non diagonal, then need to do full matvec
    for i=1:M
        B1 = func_nl(u(:,1),blank_noise,@nonlinJ, params_nl); %first nonlinear func eval
        predictor = eAdt*u(:,1) + F1*B1; %predictor

        %second nonlinear func eval and corrector
        u(:,2) = predictor + F2*(func_nl(predictor,blank_noise,@nonlinJ, params_nl)-B1);
        t(2) = t(1) + dt;
        u(:,1) = u(:,2); t(1) = t(2);
    end
end
end
