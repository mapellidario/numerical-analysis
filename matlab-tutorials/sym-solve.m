% solve equations with matlab symbolic tool
%
% define the variables
%
syms a b c real;
%
% write the equations. 
% Equations are defined as variables, in particular an equation of the form
% function(a,b,c) = 0
% becomes
% eqn(a,b,c) = function(a,b,c) - 0
%
eq1 = -b / (2*a) - 0 ;
delta = b*b - 4*a*c ;
eq2 = - delta / (4*a) + 1 / (4*a) - 2 ;
eq3 = a+b+c-5/4 ;
%
% Solve the equation
% we can use the solve() function
% this is called with
% solve(eq1,...,eqn,var1,varn)
% and returns a dictionary with the content of the variables
%
S = solve (eq1,eq2,eq3,a,b,c) ;
S.a
S.b
S.c
