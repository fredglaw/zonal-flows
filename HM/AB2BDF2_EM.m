function [u,t] = AB2BDF2_EM(init,T,M,func_noise,func_s,func_nl,params_s,params_nl)
% Implements the Adams-Bashforth Backward Differentiation Formula 2 method
% (AB2BDF3), an implicit-explicit multistep, second order accurate method.
% 
%   u_(n+1) = (4/3)u_n - (1/3)u_(n-1) + (2dt/3)A u_(n+1) +...
%             ...+ (2dt/3) [2B(u_n) - B(u_(n-1))]
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
% 
% OUTPUT: u -- numerical solution, of size N x M+1
%         t -- array of time steps, of size M+1
% 

dt = T/M; %time step

% store each step as a column vector
N = ceil(sqrt(size(init,1)));
u = zeros([size(init,1),3]); t = zeros(1,3); %initialize
u(:,1:2) = init; t(1:2) = [0,dt];
u(:,3) = u(:,2); t(3) = t(2);

A = func_s(u(:,1),params_s); %A only needs to be computed once
A = speye(size(A)) - (2*dt/3)*A; %linear system to solve at each step

% starting eval of the nonlinear function
B_prev = func_nl(u(:,1),func_noise,@nonlinJ,params_nl);

%manual choice
noise_params = params_nl(5:end);

if isdiag(A)
disp('is diag');
    B = (1./diag(A)); %SIGNIFICANT SPEED UP when A is diagonal
    for i=1:M
        B_curr = func_nl(u(:,2),func_noise,@nonlinJ,params_nl); %one nonlin func eval per step
        %form the rhs of the linear system to solve at this step
        b = (4/3)*u(:,2) - (1/3)*u(:,1) + (2*dt/3)*(2*B_curr - B_prev);
        u(:,3) = B.*b;
        
        %EM step
        curr_noise = reshape(func_noise(N,noise_params),[N*N,1]);
        u(:,3) = u(:,3) + dt*curr_noise(2:end);
        
        B_prev = B_curr; %save the nonlin func eval for next step
        t(3) = t(2) + dt;
        
        u(:,1:2) = u(:,2:3); t(1:2) = t(2:3);
    end
else
    for i=2:M
        B_curr = func_nl(u(:,2),func_noise,@nonlinJ,params_nl); %one nonlin func eval per step
        %form the rhs of the linear system to solve at this step
        b = (4/3)*u(:,2) - (1/3)*u(:,1) + (2*dt/3)*(2*B_curr - B_prev);
        u(:,3) = A\b; %solve linear system using preformed matrix
        B_prev = B_curr; %save the nonlin func eval for next step
        t(3) = t(2) + dt;
        
        u(:,1:2) = u(:,2:3); t(1:2) = t(2:3);
    end
end  
u = u(:,1:2); t = t(1:2);
end