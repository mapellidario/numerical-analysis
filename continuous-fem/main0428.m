% finite elements method for one dimensional problem
% -(cu')' = f
% u(0)=0
% c(1)u'(1)=gamma --> u'(1)=gamma/c(1)
%
clear all
close all
%
% Boundary Conditions
% Neumann non-homogeneus in x=1
gamma = 0;
%
% c and f are defined elsewhere
% define mesh
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
        N = 100;
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
    case 'random'
        % random mesh, N random points in (0,1)
        N = 20 ;
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
% contruct non uniform mesh, by transforming the current mesh
% using the point operator for term by term operations
% >> x = [1 2 3 4]
% >> x.^2;
% x: [1 4 9 16]
% x = xu.^2 ;
%

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
% last line
% Neumann omogenee
Kh(N,N-1) = -1/h(N)*c(m(N)) ;
Kh(N,N) = +1/h(N)*c(m(N)) ;
%
% >> spy(Kh)
%
% termine noto
%
fh = zeros(N,1);
for i=1:N-1
    fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
end
% utlimo elemento
fh(N) = h(N)/2*f(m(N)) + gamma;
%
% risolvo il sistema lineare
%
% uh=zeros(1:N) ;
% uh vettore colonna
uh = Kh\fh ;
% oss per al scelta delle basi, uh[i] è il valore della u(x) in x=x_i
%
%
% plot(x,uh,'o-');
% concatenate arrays: (row vector) [0 x] and (column vector) [0; uh]
% plot ([0 x], [0; uh], 'ko-')
% plot ([0 x], [0; fh], 'ok-')

% solution
figure(1)
plot([0 x], [0; uh], 'ok-');
% debug
% hold on
figure(2)
plot([0 x], [0; fh], 'ok-');

% OSS
% c=1, f=1 x=0 Dir, x=1 Neu
% In x=1 the Neumann condition is approximated! the last segment of uh
% is not perfectly horizontal
%
% Pro Tip: ctrl a, ctrl i :align the whole code!
%
% Error analysis
% (max{h})^2
% only with the code
% find a problem with a known solution, then
