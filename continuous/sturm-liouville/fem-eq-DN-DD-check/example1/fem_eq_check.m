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
clear all
% clear variables
close all


%% EXACT

syms xs real;
cs  = 10*xs ;
bs  = -xs^2  ;
as  = 1*xs+3 ;
% ues = sin(10*xs) ;
% cs  = 1+10^-10*xs ;
% bs  = 1+10^-10*xs ;
% as  = 2+10^-10*xs ;
% ues = sin(10*xs) ;
ues = sin(10*xs) + exp(-xs^2/0.01) ;
fes = -diff(cs*diff(ues,xs),xs) + bs*diff(ues,xs) + as*ues;

ue = matlabFunction(ues) ;
f = matlabFunction(fes) ;
c = matlabFunction(cs) ;
b = matlabFunction(bs) ;
a = matlabFunction(as) ;

alpha = subs(ues,xs,0) ;
beta = subs(ues, xs, 1) ;
gamma = subs(cs,xs,1)*subs(diff(ues,xs),xs,1) ;


disp(['c     : ',char(cs)])
disp(['b     : ',char(bs)])
disp(['a     : ',char(as)])
disp(['fe    : ',char(fes)])
disp(['f     : ',char(f)])
disp(['ue    : ',char(ues)])
disp(['alpha : ',char(alpha)])
disp(['beta  : ',char(beta)])
disp(['gamma : ',char(gamma)])

%
%% Boundary conditions
%
% BCtype = 'DN' ;
% u(t,0) = alpha
% c(1) u_x(t,1) = gamma
%
BCtype = 'DD' ;
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
N = 50 ;
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
fig = drawsol (uh, x, BCtype, alpha, beta) ;
hold on ;
fplot(@(x) ue(x),[0 1],'r');
title({strcat('-(c*u'')'' + b*u'' + a*u = f'),strcat('Boundary conditions: ',BCtype)}) ;
saveas(fig, strcat('fem_eq-',BCtype,'.png')) ;

%% hypothesis check

for i=1:N
    if -subs(diff(bs,xs),xs,x(i))+subs(as,xs,x(i)) <=0
        sprintf('achtung @ i: %d', i)
    end
end
