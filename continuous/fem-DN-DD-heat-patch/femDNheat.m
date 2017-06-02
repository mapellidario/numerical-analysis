%% finite elements method for one dimensional problem
%
% equation
% rho u_t -(c u_x)_x = f
%
% Curent: BC
% x = 0 : Dirichlet
% x = 1 : Neumann

% TDL: implement
% implement structures
% x = 1 : Dirichlet
% x = 0 : Neumann

%% Clean Envirnment
clear all
% clear variables
close all
%
%% Boundary conditions
%
BCtype = 'DN' ;
% u(t,0) = alpha
% c(1) u_x(t,1) = gamma
%
% BCtype = 'DD' ;
% u(t,0) = alpha
% u(t,1) = beta
%
% Boundary Conditions
% Dirichlet non-homogeneus in x=0 (homogeneus if alpha=0)
alpha = 0.1 ;
% should be ok also with alpha!=0
% Neumann   non-homogeneus in x=1 (homogeneus if gamma=0)
gamma = -1 ;
% Dirichlet non-homogeneus in x=1 ()
beta = 0.0 ;
%
%% initial condition
% initial condition
% u(0,x) = u_0 (x)
%
%% Define Mesh
%
N = 20 ;
meshtype = 'uniform' ;
% meshtype = 'random' ;
%
[x, h, m, N] = makemesh(meshtype, N) ;
%
%% Matrices
%
% This part depends on the BC implicitly, meaning it "only" changes the
% stricture of the matrices
%
% diffusion matrix
%
[Kh] = makeK (N, h, m, BCtype) ;
%
% reaction matirx
%
[Mh] = makeM (N, h, m, BCtype) ;
%
%% Constant Term
%
% this part dependends on BC explicitly, both Dirichlet and Neumann!
[fh] = makef (N, h, m, BCtype, alpha, beta, gamma) ;
%
%
%% Solving ODE
% Solve ODE system with implicit Euler method
% From here on we have the ODE system with Kh Mh and fh
%
% we chose dt, setting total time to analyse and the number of desired
% steps, so that we can coherently analyse the evolution changing the
% numbers of time iterations
%
% final time and number of time steps (kmax)
time = 1 ;
kmax = 20 ;
%
dt = time / kmax ;
%
% this part is BC independent!
[uh, uh0] = solveode (N, x, dt, kmax, Kh, Mh, fh, BCtype) ;
%
%% Draw Solution
% plot
%
% this depends on Dirichlet BC!
drawsol   (uh, uh0, x,        kmax, BCtype, alpha, beta) ;
drawpatch (uh, uh0, x, N, dt, kmax, BCtype, alpha, beta) ;
