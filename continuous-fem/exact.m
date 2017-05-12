clear all
close all

% write a problem with a known solution
%
%
syms x real;
%
% chose c: c(x)=1+x^2
%
c = 1+x^2;
%
% chose exact solution
% with the only constrain that ue(0)=0
ue = sin(5*x)+log(1+x) ;
%
%
% -(cu')' = f
% c has to be regular for having convergence of order h^2
%
% given c and given u, we find f (with symbolic calculus)
% we now have c,u,f such that we know an exact solution!
%
% compute f and gamma
% f= -(cu')'
% gamma = c(1)ue'(1)
%
% >> syms x real
% >> c = 1+x^2
% >> diff(c,x)
% >> int(c,x,0,1)

uep = diff(ue,x) ;
f = -diff(c*uep,x) ;
% or also: f = -diff(c*diff(ue,x),x);
% NB this can be done if c and f are regular! ()
%
% how to evaluate a funcion in a point?
% subs() command to change variable
% >> subs(c,x,x+1)
% >> subs(c,x,1)
% but this is still symbolic
% eval() evaluate the symbolic expression and returns a real value
% >> eval(subs(c,x,1))

gamma = subs(c,x,1)*subs(uep,x,1) ;

% now we set c, f ang gamma to the result we obtained here
% so that we have a method to compare uh and ue!

%
% >> ezplot('sin(5*x)+log(1+x)',[0 1])
