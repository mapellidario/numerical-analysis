clear all
close all
%
% define variables for two generic parabolas:
%
syms x A B C D E F real;
syms xs c1 c2 f1 f2 real;
%
u1 = A*x^2 + B*x + C ;
u2 = D*x^2 + E*x + F ;
%
% equation in x = 0.
% Achtung: equation are written in f(x) = 0 form
% (=0 is hidden when we write eqn)
%
eq1 = subs(u1,x,0) - 0 ;
%
% differential equation in ]0,xs[
%
eq2 = -diff(c1*diff(u1,x),x) - f1 ;
%
% continuity of u in xs
%
eq3 = subs(u1,x,xs) - subs(u2,x,xs);
%
% continuity of w=c*u' in xs
%
% stress in zone 1
w1 = c1 * diff(u1,x);
% stress in zone 2
w2 = c2 * diff(u2,x);
%
eq4 = subs(w1,x,xs) - subs(w2,x,xs);
%
% differential equation in ]xs,1[
%
eq5 = -diff(c2*diff(u2,x),x) - f2 ;
%
% equation in 1
%
% Dirichlet
eq6 = subs(u2,x,1) - 0;
% Neumann
% eq6 = subs(c2*diff(u2,x),x,1)-0;

%
% Solve system of equations
% solve($equations, $variables) ;
%
S = solve(eq1, eq2, eq3, eq4, eq5, eq6, A, B, C, D, E, F) ;
%
% check the solutions
% >> S.A
% >> S.B
%
% Plot the solution
% Eliminate S. so that the parameters name is nicer.
% This is accomplished by assigning the value of the slution contained
% into the structure S.to simple variables, redefining the symbolic
% variables to normal variables
%
A = S.A;
B = S.B;
C = S.C;
D = S.D;
E = S.E;
F = S.F;
%
% u1 and u2 also nees to be redefined in terms of the new
% simple variebles and not of the symbolic real variables
%
u1 = A*x^2+B*x+C;
u2 = D*x^2+E*x+F;
%
% define some handy parameters that will be used to control some
% characteristics of the plot
%
% first set of interesting parameters
vc1 = 1;
vc2 = 1;
vf1 = 1;
vf2 = 100;
vxs = 0.5;

% % second set of interesting parameters values
% vc1 = 1;
% vc2 = 2;
% vf1 = 1;
% vf2 = 1;
% vxs = 0.5;
%
% this set is used to check that the solution makes sense:
% if c and f are continuous (vc1=vc2 and vf1=vf2), then should happen
% that u1 = u2 (A=D, B=E, C=F), and it does not depend by xs
% vc1 = 1;
% vc2 = 2;
% vf1 = 1;
% vf2 = 1;
% vxs = 0.5;
%
%
u1 = subs(u1,[c1 c2 f1 f2 xs], [vc1 vc2 vf1 vf2 vxs]);
u2 = subs(u2,[c1 c2 f1 f2 xs], [vc1 vc2 vf1 vf2 vxs]);
%
% plot u1 and u2 on the same plot
% create figure
figure(1)
fplot(u1,[0 vxs])
% force persisten old graph
hold on
% new graph, on the same figure of the old
fplot(u2,[vxs 1])
