% finite elements method
clear all
close all
%
% c and f are defined elsewhere
% ------------ DEFINE MESH ---------------
%
% we start with a uniform mesh
% >> linspace(xmin, xmax, M)
% divide the interval [xmin,xmax] in M-1 intervals
%
% we want N intervals
%
N = 10;
x = linspace(0,1,N+1);
%
% problem: x_1=0, x_N+1 = 1
%
x = x(2:end);
%
% this new x is the old x, but starts from 2nd element
% this aligns matlab vector with the one used during lectures
%
% compute h and m
%
h = zeros(1,N);
m = zeros(1,N);
%
% i-1 = 0 for i=1, so this array starts
%
for i=2:N
        h(i) = x(i)-x(i-1);
        m(i) = (x(i-1)+x(i))/2 ;
end
%
% i=1
h(1) = x(1)-0;
m(1) = x(1)/2;
%
% -------------- MATRIX KH ------------
%
Kh = zeros(N,N);
%
for i=1:N-1
  %
  if i>1
    % this element is not present in the first row
    Kh(i,i-1) = -1/h(i)*c(m(i));
  end
  %
  Kh(i,i) = 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1));
  %
  Kh(i,i+1) = -1/h(i+1)*c(m(i+1));
  %
end
%
% last row
Kh(N,N-1) = -1/h(N)*c(m(N));
Kh(N,N) = 1/h(N)*c(m(N));
%
% run the program and
% >> whos
% >> spy(Kh)
%
%
% --------- CONSTANT TERM --------------
%
fh = zeros(N,1);
%
for i=1:N-1
  fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1));
end
%
% last element
%
fh(N) = h(N)/2*f(m(N));
%
%
% ------------ SOLVE LINEAR SYSTEM ------------
%
uh = Kh\fh;
% because of choice of basis (Lagrange)
% the value of uh in a point is the value of the function in that point
%
% plot: 'k-' graphic option
% be careful when concatenating row and column vectors (use ; if column)
plot([0 x], [0; uh], 'k-');
