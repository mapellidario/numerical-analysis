%
% matlab symbolic calculus: generic operations
% inside symbolic toolbox
%
syms x real
expand ((x+3)^3)
int (x^2, 0, 2)


%
% define a generic parabola:
%
syms x A B C real;
u1 = A*x^2 + B*x + C ;

% >> u1
% substitute to a simbolic variable another variable or a defined value
% >> subs(u1,x,0)
% differentiate
% >> diff(u1.x)
% integrate
% >> int(u1,x)
% >> int(u1, x,0,1)
% >> int (u1^3,x,0,1)
