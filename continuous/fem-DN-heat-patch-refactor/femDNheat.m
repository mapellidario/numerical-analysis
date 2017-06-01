% finite elements method for one dimensional problem
%
% equation
% rho u_t -(c u_x)_x = f
%
% Boundary conditions
% u(t,0) = alpha
% c(1) u_x(t,1) = gamma
%
% initial condition
% u(0,x) = u_0 (x)
%
clear all
% clear variables
close all
%
% Boundary Conditions
% Dirichlet non-homogeneus in x=0 (homogeneus if alpha=0)
alpha = 0 ;
% should be ok also with alpha!=0
% Neumann   non-homogeneus in x=1 (homogeneus if gamma=0)
gamma = 0 ;

%
% c and f are defined elsewhere
%% ------------ DEFINE MESH ---------------
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
        N = 100 ;
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
        N = 50 ;
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

%
% Heat modificatons starts from here
%


%% -------------- MATRIX KH ------------
% diffusion matrix
%
% Kh = zeros(N,N);
Kh = sparse(N,N);
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
% reaction matirx
%
% Mh = zeros(N,N);
Mh = sparse(N,N);
for i=1:N-1
    %
    if i>1
        % only from second row
        Mh(i,i-1) = - h(i)/4*rho(m(i)) ;
    end
    %
    Mh(i,i) = + h(i)/4*rho(m(i)) + h(i+1)/4*rho(m(i+1)) ;
    %
    Mh(i,i+1) = h(i+1)/4*rho(m(i+1)) ;
end
% last row
% Neumann homogenous
Mh(N,N-1) = h(N)/4*rho(m(N)) ;
Mh(N,N)   = h(N)/4*rho(m(N)) ;
%
% >> whos
% >> spy(Mh)
%
%% --------- CONSTANT TERM --------------
%
fh = zeros(N,1);
for i=1:N-1
    fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
end
% first element for dirichlet in x=0
fh(1) = fh(1) + alpha/h(1)*c(m(1)) - alpha*h(1)/4*rho(m(1));
% last element for neumann in x=1
fh(N) = h(N)/2*f(m(N)) + gamma;

%
%% ---- TIME ------
% Solve ODE system with implicit Euler method
% From here on we have the ODE system with Kh Mh and fh
%
% we chose dt, setting total time to analyse and the number of desired
% steps, so that we can coherently analyse the evolution changing the
% numbers of time iterations
%
% final time and number of time steps (kmax)
time = 1 ;
kmax = 100 ;
%
dt = time / kmax ;
%
% we now use implicit euler. we save every step
uh = zeros(N,kmax) ;
% uh(;,k) is the solution at step k.
%
% OSS: we do not have the initial condition in uh(:,0) !
% we have to define the initial condition in a separate column vector
uh0 = zeros(N,1) ;
for i=1:N
    % interpolating initial condition u0
    % uh0(i) = u0(x(i)) ;
    uh0(i) = rand() ;
end
%
% horrible, but we have to pay that vectors start from 1 in matlab
% and we want to adhere to conventions used during lectures
% implicit euler step
for k=0:kmax-1
    % uh0 = uh(:,0) -> uh(:,1)
    if k==0
        uh(:,k+1) = (1/dt*Mh + Kh)\(1/dt*Mh*uh0+fh) ;
    else
        % uh0 = uh(:,k) -> uh(:,k+1)
        % step of implicit Euler
        uh(:,k+1) = (1/dt*Mh + Kh)\(1/dt*Mh*uh(:,k)+fh) ;
    end
end
%
%% ----- DRAW SOLUTION -------
% plot
figure(1)
%
% initial condition
plot ([0 x], [alpha;uh0], 'ro-' )
%
% plot the solution in every time step
hold on
for k=1:kmax
    plot ([0 x], [alpha; uh(:,k)], 'b' )
end


figure(2)
% draw the solution as a surface
% we use the patch command, which is supergeneric!
% being the grid cartesian we can use a more specific command, but
% this command is pretty generic and is nice to use it

% we assume that there is Neumann in x=1,
for i=1:N-1
    for k=1:kmax-1
        % point 1
        t1 = k*dt ;
        x1 = x(i) ;
        u1 = uh(i,k) ;
        % point 2
        t2 = (k+1)*dt ;
        x2 = x1 ;
        u2 = uh(i,k+1) ;
        % point3
        t3 = t2 ;
        x3 = x(i+1) ;
        u3 = uh(i+1,k+1) ;
        % point 4
        t4 = t1 ;
        x4 = x3 ;
        u4 = uh(i+1,k) ;
        %
        % uh in scala di colore
        % patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4]) ;
        % la voglio tutta bianca, quindi le varie u le considera come
        % altezza
        patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;
        
    end
end

% t=0
for i=1:N-1
    % point 1
    t1 = 0 ;
    x1 = x(i) ;
    u1 = uh0(i) ;
    % point 2
    t2 = dt ;
    x2 = x1 ;
    u2 = uh(i,1) ;
    % point3
    t3 = t2 ;
    x3 = x(i+1) ;
    u3 = uh(i+1,1) ;
    % point 4
    t4 = t1 ;
    x4 = x3 ;
    u4 = uh0(i+1) ;
    %
    % uh in scala di colore
    % patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4]) ;
    % la voglio tutta bianca, quindi le varie u le considera come
    % altezza
    patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;
end


% x=0
for k=1:kmax-1
    % point 1
    t1 = k*dt ;
    x1 = 0 ;
    u1 = alpha ;
    % point 2
    t2 = (k+1)*dt ;
    x2 = x1 ;
    u2 = alpha ;
    % point3
    t3 = t2 ;
    x3 = x(1) ;
    u3 = uh(1,k+1) ;
    % point 4
    t4 = t1 ;
    x4 = x3 ;
    u4 = uh(1,k) ;
    %
    % uh in scala di colore
    % patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4]) ;
    % la voglio tutta bianca, quindi le varie u le considera come
    % altezza
    patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;
    
end

% x=0, t=0
% point 1
t1 = 0 ;
x1 = 0 ;
u1 = alpha ;
% point 2
t2 = dt ;
x2 = x1 ;
u2 = alpha ;
% point3
t3 = t2 ;
x3 = x(1) ;
u3 = uh(1,1) ;
% point 4
t4 = t1 ;
x4 = x3 ;
u4 = uh0(1) ;
patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;
