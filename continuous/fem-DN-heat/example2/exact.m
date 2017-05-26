% Problem -(cu')' + au = f
% BC Dirichlet in x=0
% BC Neumann in x=gamma
% Ater fixing a solution, we compute the costant term,
% given c and f
% The purpose is having an exact solution to evaluate the
% error of the fem
%
clear all
close all
%
syms x real;
%
% choose c:
% c has to be regular for having convergence of order hmax^2
% c has to be c(x) >= c_0 > 0
% this requires that in
% c.m:
% function y=c(x)
% y=1+x^2
c = 1+x^2 ;
%
% choose ca:
% a has to be regular and non-negative
% a has to be a(x) >= 0
% this requires that in
% a.m:
% function y=a(x)
% y=1+x^2
a = (sin(x))^2 ;
%
% choose exact solution
% this solution requires
% femDN.m :
% alpha=0
% this solution has to be written into
% ue.m
%
ue = cos(5*x)+log(1+x) ;
%
% given c and given u, we find f (with symbolic calculus)
% we now have c,u,f such that we know an exact solution!
%
% compute f and gamma
% f= -(cu')'
% gamma = c(1)ue'(1)
%
uep = diff(ue,x) ;
fe = -diff(c*uep,x) + a*ue;
% or also in one shot: f = -diff(c*diff(ue,x),x) + a*ue;
% NB this can be done if c and ue are regular! ()
%
% how to evaluate a funcion in a point?
% subs() command to change variable
% >> subs(c,x,x+1)
% >> subs(c,x,1)
% but this is still symbolic
% eval() evaluate the symbolic expression and returns a real value
% >> eval(subs(c,x,1))

alpha = subs(ue,x,0) ;
gamma = subs(c,x,1)*subs(uep,x,1);

disp(['c     (c.m)    : ',char(c)])
disp(['a     (a.m)    : ',char(a)])
disp(['fe    (f.m)    : ',char(fe)])
disp(['ue    (ue.m)   :  ',char(ue)])
disp(['alpha (femDN.m): ',char(alpha)])
disp(['gamma (femDN.m): ',char(gamma)])

% now we set c, f, alpha and gamma to the result we obtained here
% so that we have a method to compare uh and ue!
