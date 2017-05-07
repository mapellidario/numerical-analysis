% finite elements method
clear all
close all
%
% c and f are defined elsewhere
% define mesh
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
% i-1 = 0 for i=1, so this cycle starts from 2
%



for i=2:N
        h(i) = x(i)-x(i-1);
        m(i) = (x(i-1)+x(i))/2 ;
end
%
Kh = zeros(N,N);
%
for i=1:N-1
    %
    if i>1
        % first line
        Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
    end
    %
    Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
    %
    Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
end
% last line
Kh(N,N-1) = -1/h(N)*c(m(N)) ;
Kh(N,N) = +1/h(N)*c(m(N)) ;

