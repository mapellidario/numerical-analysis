function [hmax,errmax]=fem_eq_func(N, meshtype, c, b, a, f, ue, BCtype, alpha, beta, gamma, bs, as, xs) ;

disp(['N ', num2str(N)])

%% finite elements method for one dimensional problem
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
% * This is similar to fem_eq.m, but it is target to check the fem
% method. This automatically set the correct parameters to compute
% the approximate solution, given the exact one.
% We can check the validity of the code by setting c,b,a and an exact
% solution, and then compute f.
% Then we use f to compute the approximate solution and
% plot it together with the exact solution.
% This process is here done automatically
% Note that it should work both with DD and DN BCs.
%
% * hypotesis check
%
% fem works if $ -b_x/2+a > 0 $, and atthe end of the program we check
% if this condition holds. If it does not, we can not trust the fem
% to properly work!


%% Clean Envirnment
% clear all
% clear variables
% close all

%
%% Boundary conditions
%
% BCtype = 'DN' ;
% u(t,0) = alpha
% c(1) u_x(t,1) = gamma
%
% BCtype = 'DD' ;
% u(t,0) = alpha
% u(t,1) = beta
%
% Boundary Conditions
% Dirichlet non-homogeneus in x=0 (homogeneus if alpha=0)
% alpha = 0 ;
% should be ok also with alpha!=0
% Neumann   non-homogeneus in x=1 (homogeneus if gamma=0)
% gamma = 1 - 10*sin(5) ;
% gamma = 0 ;
% Dirichlet non-homogeneus in x=1 ()
% beta = cos(5) + log(2) ;
% beta = 0 ;
%
%% initial condition
% initial condition
% u(0,x) = u_0 (x)
%
%% Define Mesh
%
% N = 1000 ;
% meshtype = 'uniform' ;
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
[Kh] = makeK (N, h, m, c, BCtype) ;
%
% transport term
%
[Bh] = makeB (N, h, m, b, BCtype) ;
%
% reaction matirx
%
[Ah] = makeA (N, h, m, a, BCtype) ;
%
%% Constant Term
%
% this part dependends on BC explicitly, both Dirichlet and Neumann!
[fh] = makef (N, h, m, BCtype, c,b,a,f, alpha, beta, gamma) ;
%
%

%% solve linear system

uh = (Kh + Bh + Ah) \ fh ;

%% Draw Solution
% plot
%
% this depends on Dirichlet BC!
% drawsol (uh, x, BCtype, alpha, beta) ;

% hold on ;
% fplot(@(x) ue(x),[0 1],'r');

%% hypothesis check

for i=1:N
    if -subs(diff(bs,xs),xs,x(i))+subs(as,xs,x(i)) <=0
        sprintf('achtung @ i: %d', i)
    end
end

%% return

errmax=0;
switch strcat(BCtype)
    case 'DN'
        for i=1:N
            erri = abs(uh(i)-ue(x(i))) ;
            if erri>errmax
                errmax=erri;
            end
        end
    case 'DD'
        for i=1:N-1
            erri = abs(uh(i)-ue(x(i))) ;
            if erri>errmax
                errmax=erri;
            end
        end
end        
        % get maximum of the array h
hmax = max(h) ;

end
