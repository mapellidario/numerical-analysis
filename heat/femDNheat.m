% finite elements method for one dimensional problem
% rho u_t - (cu')' = f
%
% boundary conditions
% u(t,0)=alpha
% c(1)u_x(t,1)=gamma
%
% initial condition
% u(0,x) = u0(x)
%
clear all
close all
%
% Boundary Conditions
% Dirichlet non-homogeneus in x=0 (homogeneus if alpha=0)
alpha = 0 ; 
% Neumann non-homogeneus in x=1 (homogeneus if gamma=0)
gamma = 0 ;

%
% c and f are defined elsewhere
% ------------ DEFINE MESH ---------------
%
 mesh = 'uniform' ;
% mesh = 'random' ;
%
switch mesh
    case 'uniform'
        % we start with a uniform mesh (x_uniform)
        % >> linspace(xmin, xmax, M)
        % divide the interval [xmin,xmax] in M-1 intervals (with M points)
        %
        % we want N intervals
        N = 10;
        xu = linspace(0,1,N+1);
        %
        % problem: x_1=0, x_N+1 = 1
        %
        xu = xu(2:end);
        %
        x = xu ;
        %
        % this new x is the old x, but starts from 2nd element
        % this aligns matlab vector with the one used during lectures
        % in matlab vectors strart from 1, there is no index = 0 in arrays.
    case 'transformed'
        % contruct non uniform mesh, by transforming the current mesh
        % using the point operator for term by term operations
        % >> x = [1 2 3 4]
        % >> x.^2;
        % x: [1 4 9 16]
        % x = xu.^2 ;
        %
        disp('not implemented yet')
        %
    case 'random'
        % random mesh, N random points in (0,1)
        N = 200 ;
        %
        x = rand(1,N) ;
        % sort the array, from lower to higher
        x = sort(x) ;
        % add x=0 and x=1
        x = [0 x 1] ;
        % remove repeated elements
        x = unique(x) ;
        %
        % remove x=0
        x = x(2:end) ;
        % redefine N
        N = length(x) ;
    otherwise
        %
        disp('wrong mesh name')
        %
end
%
% compute h and m
% this works for a generic mesh
%
h = zeros(1,N);
m = zeros(1,N);

%
% i-1 = 0 for i=1, so this cycle starts from 2
%
h(1) = x(1) - 0;
m(1) = (x(1) + 0)/2;
for i=2:N
        h(i) = x(i)-x(i-1);
        m(i) = (x(i-1)+x(i))/2 ;
end
%
% inspect the mesh
% >> bar(h)
%
% -------------- MATRIX KH ------------
%
Kh = zeros(N,N);
%
for i=1:N-1
    %
    if i>1
        % only from second row
        Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
    end
    %
    Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
    %
    Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
end
% last row
% Neumann homogenous
Kh(N,N-1) = -1/h(N)*c(m(N)) ;
Kh(N,N) = +1/h(N)*c(m(N)) ;
%
%
% -------------- MATRIX MH ------------
%
Mh = zeros(N,N);
%
for i=1:N-1
    %
    if i>1
        % only from second row
        Mh(i,i-1) = - h(i)/4*rho(m(i)) ;
    end
    %
    Mh(i,i) = h(i)/4*rho(m(i)) + h(i+1)/4*rho(m(i+1)) ;
    %
    Mh(i,i+1) = h(i+1)/4*rho(m(i+1)) ;
end
% last row
% Neumann homogenous
Mh(N,N-1) = h(N)/4*rho(m(N)) ;
Mh(N,N) = h(N)/4*rho(m(N)) ;
%
%
% >> whos
% >> spy(Kh)
% >> spy(Mh)
%
% --------- CONSTANT TERM --------------
%
fh = zeros(N,1);
for i=1:N-1
    fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
end
% first element for dirichlet in x=0
fh(1) = fh(1) + alpha/h(1)*c(m(1)) - alpha*h(1)*rho(m(1));
% last element for neumann in x=1
fh(N) = h(N)/2*f(m(N)) + gamma;
%
% ------------ SOLVE ODE SYSTEM with EXPLICIT EULER METHOD ------------
% ODE with Mh, Kh, fh
%
% choose dt
dt = 0.1;
%
% choose kmax (n of steps in time)
kmax = 100;
%%
uh = zeros(N, kmax);
%
% uh(:,k) is the solution at the k-th step
%
% define initial condition uh0 (column vector with initial values)
%
uh0 = zeros(N,1);
%
% interpolation of initial condition u0
%
for i=1:N
    uh0(i) = u0(x(i));
end
%
% explicit Euler step
%
for k=0:kmax-1
    if k==0 % uh(:,0) is not allowed in matlab so i use uh0
        uh(:,k+1) = (1/dt*Mh + Kh)\(1/dt*Mh*uh0+fh);
    else % uh(:,k) exists
        uh(:,k+1) = (1/dt*Mh + Kh)\(1/dt*Mh*uh(:,k)+fh);
    end
end
%
%
%
% -------------- DRAW UH  --------------
figure(1);
%
% draw initial condition
%
plot([0 x],[alpha; uh0],'r') 
%
% draw uh in all the time steps
%
hold on
for k=1:kmax
    plot([0 x],[alpha; uh(:,k)],'b') 
end