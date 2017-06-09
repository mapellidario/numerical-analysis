% Relation between hmax and error in the norm of the maximum
% This program is the same as the one used for fem_eq_exact.m , but it is
% not used to validate the fem, but to estimate the error of the method.
%
% * equation
%
% -(c*u_x)_x + b*u_x + a*u = f
% -(c*u')' + b*u' + a*u = f
%
% * BC
%
% DN
% x = 0 : Dirichlet
% x = 1 : Neumann
% u(0) = alpha
% c(1)*u_x(1) = gamma
%
% DD
% x = 0 : Dirichlet
% x = 1 : Dirichlet
% u(0) = alpha
% u(1) = beta
%
% * parameters
%
% c: set in c.m
% b: set in b.m
% a: set in a.m
% f: set in f.m
%
% * mesh
%
% Set mesh characteristics: meshtype and number of mesh points
% 
% 
%     Linear model Poly1:
%     f(x) = p1*x + p2
%     Coefficients (with 95% confidence bounds):
%       p1 =       2.001  (2, 2.001)
%       p2 =       1.805  (1.804, 1.806)
% We note that there is a quadratic relation, since p1 is pretty close to 2
%
clear all
close all

%% EXACT

syms xs real;
cs  = 10*xs ;
bs  = -xs^2 ;
as  = 1*xs + 3 ;
% cs  = 1+10^-10*xs ;
% bs  = 1+1^-10*xs ;
% as  = 1+10^-10*xs ;
% ues = exp(-xs^2/0.01) ;
ues = sin(10*xs) + exp(-xs^2/0.01) ;
fes = -diff(cs*diff(ues,xs),xs) + bs*diff(ues,xs) + as*ues;

ue = matlabFunction(ues) ;
f = matlabFunction(fes) ;
c = matlabFunction(cs) ;
b = matlabFunction(bs) ;
a = matlabFunction(as) ;

alpha = subs(ues,xs,0) ;
beta = subs(ues, xs, 1) ;
gamma = subs(c,xs,1)*subs(diff(ues,xs),xs,1) ;

disp(['c     : ',char(cs)])
disp(['b     : ',char(bs)])
disp(['a     : ',char(as)])
disp(['fe    : ',char(fes)])
disp(['f     : ',char(f)])
disp(['ue    : ',char(ues)])
disp(['alpha : ',char(alpha)])
disp(['beta  : ',char(beta)])
disp(['gamma : ',char(gamma)])

%% cicle

N = 20:10:100 ;
lN = length(N) ;

meshtype='uniform' ;
% meshtype='random' ;

BCtype='DN' ;
% BCtype='DD' ;

h=zeros(lN ,1);
err=zeros(lN ,1);

for i = 1:lN
    [h(i),err(i)]=fem_eq_func(N(i), meshtype, c, b, a, f, ue, BCtype, alpha, beta, gamma, bs, as, xs) ;
end

%% plot

% disp(['c     : ',char(cs)])
% disp(['b     : ',char(bs)])
% disp(['a     : ',char(as)])
% disp(['fe    : ',char(fes)])
% disp(['f     : ',char(f)])
% disp(['ue    : ',char(ues)])
% disp(['alpha : ',char(alpha)])
% disp(['beta  : ',char(beta)])
% disp(['gamma : ',char(gamma)])

FigLogLog = figure(1);
loglog(h,err,'o-');
grid on
switch strcat(BCtype)
    case 'DN'
        title({'-(cu'')''+bu''+au=f';'Dirichlet - Neumann'}); 
    case 'DD'
        title({'-(cu'')''+bu''+au=f';'Dirichlet - Dirichlet'}); 
end
xlabel('err');
ylabel('hmax');

saveas(FigLogLog, strcat('fem_eq_error-',BCtype,'.png')) ;

logh = log(h) ;
logerr = log(err) ; 
f=fit(logh,logerr,'poly1') ; 
FitLogLog = figure(2);
switch strcat(BCtype)
    case 'DN'
        title({'-(cu'')''+bu''+au=f - loglogfit';'DN'})
    case 'DD'
        title({'-(cu'')''+bu''+au=f - loglogfit';'DD'})
end
plot(f,logh,logerr,'o') ;
switch strcat(BCtype)
    case 'DN'
        title({'-(cu'')''+bu''+au=f';'Dirichlet - Neumann'}) ;  
    case 'DD'
        title({'-(cu'')''+bu''+au=f';'Dirichlet - Dirichlet'}) ;  
end
ylabel('log(err)') ;
xlabel('log(hmax)') ;
saveas(FitLogLog, strcat('fem_eq_error1-',BCtype,'.png')) ;

disp(['Nmin ',num2str(N(1)), ' Nmax ', num2str(N(lN))])