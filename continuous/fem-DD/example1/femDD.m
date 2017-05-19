% finite elements method for one dimensional problem
% -(cu')' = f
% u(0)=alpha
% u(1)=beta
%
clear all
close all
%
% Boundary Conditions
% Dirichlet non-homogeneus in x=0
alpha = 0 ;
% Dirichlet non-homogeneus in x=1
beta = 0 ;
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
Kh = zeros(N-1,N-1);
%
for i=1:N-1
    %
    if i>1
        % this term is not present in the first row
        Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
    end
    %
    Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
    %
    if i<N-1
        % this term is not present in the last row. meditate on this!
        Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
    end
end
%
% >> spy(Kh)
%
% termine noto, vettore colonna
%
fh = zeros(N-1,1);
for i=1:N-1
    fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
end
% prima riga per dirichlet non omogeneo in x=0
fh(1) = fh(1) + alpha/h(1)*c(m(1)) ;
% ultima riga per dirichlet non omogeneo in x=1
fh(N-1) = fh(N-1) + beta/h(N)*c(m(N)) ;

%
% risolvo il sistema lineare
%
% uh=zeros(1:N-1) ;
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
plot([0 x], [alpha; uh; beta], 'ok-');
% scale the axis
% axis([xmin xmax ymin ymax])
% axis([0 1 -1 4]) ;


% debug
% hold on
% figure(2)
% plot([0 x], [0; fh], 'ok-');

% % compare to exact solution
% hold on ;
% %fplot('sin(5*x)+log(1+x)',[0 1],'r');
% fplot(@(x)sin(5*x)+log(1+x),[0 1],'r');
%
% % OSS points of uh are not on u! it is not an interpolating solution!
% % in can be shown that if h->0, then the distance from the point of uh
% % and the correponing points of u gets lower.
% % we want to estimate how quickly the method converges.
%
% % there are various methods to estimate the error
% % Here we consider the maximum error at the nodes
% errmax=0;
% for i=1:N
%     erri = abs(uh(i)-ue(x(i))) ;
%     if erri>errmax
%         errmax=erri;
%     end
% end
% % get maximum of the array h
% hmax = max(h) ;
% format shorte ;
% [N hmax errmax]

% OSS: asintotic behaviour of the error
% uniform mesh
% N hmax errmax
%    1.0000e+01   1.0000e-01   2.1264e-02
%    2.0000e+01   5.0000e-02   5.2103e-03
%    1.0000e+02   1.0000e-02   2.0708e-04
% random mesh
% the error is dominated by the larger h!
% N hmax errmax
%    2.1000e+01   1.5985e-01   7.8718e-02
%    2.0100e+02   3.4237e-02   6.5123e-04

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
