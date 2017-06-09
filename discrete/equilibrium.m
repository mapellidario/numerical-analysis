%% Discrete model of elasticity.
% system of masses and spring hanging from the ceiling

% erase all previous variables
clear all
% close all previous figures and open files
close all
%
%% System Parameters
%
% number of masses
N = 300 ; 
%
% number of springs
endcondition = 'free' ;
% endcondition = 'spring' ;

switch strcat(endcondition)
    case 'free'
        % free end 
        M = N ;
    case 'spring'
        % spring at the end
        M= N + 1 ;
end 

%% first create the matrix, giving it a proper size
A = zeros(M,N);
%
% fill in the diagonal
for i=1:N
  A(i,i) = 1;
end
%
% fill in the subdiagonal
for i=2:M
  A(i,i-1) = -1;
end
%
% check that the matrix A is correct
% >> l_03_17_0
% >> A

%% define matrix C
%
C = eye(M,M);
% fill in the last half of the diagonal with a different c;
for i=round(M/2):M
  C(i,i) = 10;
end
%
% check that the matrix C is correct
% >> l_03_17_0
% >> C
%
%% create matrix K
%
K = A' * C * A;
% >> check matrix K
% >> pretty(K)
%
% syms c1 c2 c3 c4 c5 c6 real
% C = diag ([c1 c2 c3 c4 c5 c6]);
%
%% create force vector
%
f = ones(N,1);
%
%% solve the linear system
%
u = K\f;

%% draw the solution
%
graph = figure(1) ;
bar(u);

% save the image
title({'discrete model', strcat(endcondition)}) ;
imagename='discrete' ;
saveas(graph,strcat(imagename,'-', endcondition,'.png'));
% 
% notice that the result are parabola arcs, connected C^0
