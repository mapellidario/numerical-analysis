%% Finite elements method for one dimensional problem
%
% * equation
%
% -(c*u_x)_x + b*u_x + a*u = f
% -(c*u')' + b*u' + a*u = f
%
% * BC
%
% DN
% x = 0 : Dirichlet
% x = 1 : Neumann
% u(0) = alpha
% c(1)*u_x(1) = gamma
% 
% DD
% x = 0 : Dirichlet
% x = 1 : Dirichlet
% u(0) = alpha
% u(1) = beta
%
% * parameters
%
% c: set in c.m
% b: set in b.m
% a: set in a.m
% f: set in f.m
%
% * mesh
%
% Set mesh characteristics: meshtype and number of mesh points
%
% * Exact solution
%
% We can check the validity of the code by setting c,b,a and an exact
% solution, and then compute f.
% Then we use f to compute the approximate solution and 
% plot it together with the exact solution
%
% * TDL
%
% 1. implement structures
% 2. x = 1 : Dirichlet
% 3. x = 0 : Neumann

%% Clean Envirnment
%

clear all
% clear variables
close all

%% Boundary conditions
%

BCtype = 'DN' ;
% BCtype = 'DD' ;

% Dirichlet non-homogeneus in x=0 (homogeneus if alpha=0)
alpha = 0 ;
% Neumann   non-homogeneus in x=1 (homogeneus if gamma=0)
gamma = 100*cos(10) ;
% Dirichlet non-homogeneus in x=1 ()
beta = sin(10) ;

%% Define Mesh
%

N = 100 ;
%
meshtype = 'uniform' ;
% meshtype = 'random' ;
%
[x, h, m, N] = makemesh(meshtype, N) ;

%% Matrices
%
% Depending on the BCtype, this creates the correct matrices

% diffusion matrix
%
[Kh] = makeK (N, h, m, BCtype) ;
%
% transport term
%
[Bh] = makeB (N, h, m, BCtype) ;
%
% reaction matirx
%
[Ah] = makeA (N, h, m, BCtype) ;

%% Constant Term
%
% This creates the right constant term

[fh] = makef (N, h, m, BCtype, alpha, beta, gamma) ;

%% solve linear system
%

uh = (Kh + Bh + Ah) \ fh ;

%% Draw Solutionc
% 
% Plot the approximated solution.
% Be careful, uh is not defined where Dirichlet BC are set!
% We have to set them values manually inside the plot
%
% In order to check the validity of the method, we can also plot the
% exact solution

fig = drawsol (uh, x, BCtype, alpha, beta) ;
hold on ;
fplot(@(x) ue(x),[0 1],'r');
title({'-(c*u'')'' + b*u'' +a*u = f',strcat('Boundary conditions: ',BCtype)}) ;
saveas(fig, strcat('fem_eq-',BCtype,'.png')) ;