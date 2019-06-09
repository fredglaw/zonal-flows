close all; clear all;
% Sample script to test the accuracy on nonlinJ, on periodic functions
N = 64; %number of nodes in each direction
L = 4*pi; %width of interval
sc = L / (2*pi); %scaling factor
x = linspace(-L/2,L/2,N+1); x(end) = [];

[X,Y] = meshgrid(x,x); %build square grid

%test functions
f = sin(X).*sin(Y);
g = sin(X).*sin(2*Y);

%compute (f_x)(g_y)-(f_y)(g_x)
truth = (-cos(X).*sin(Y)).*(-2*cos(2*Y).*sin(X)) - (-sin(X).*cos(Y)).*(-cos(X).*sin(2*Y));
temp = ifft2(nonlinJ(fft2(f), fft2(g),sc));
disp(norm(imag(temp)));
temp_r = real(temp);

disp(['max deviation from truth is ', num2str(norm(temp_r-truth,Inf))]);

figure(1);
surf(truth);
title('truth');

figure(2);
surf(temp_r);
title('numerical soln');