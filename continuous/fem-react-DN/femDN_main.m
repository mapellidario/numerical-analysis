%
% Relation between hmax and error in the norm of the maximum
% This program uses femDN_func.m, which is a modified version of femDN.m
%
% we consider finite elements method for one dimensional problem
% -(cu')' + au = f
% u(0)=alpha
% c(1)u'(1)=gamma --> u'(1)=gamma/c(1)
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

N = 10:10:500 ;
lN = length(N) ;

h=zeros(lN ,1);
err=zeros(lN ,1);

for i = 1:lN
    [h(i),err(i)]=femDN_func(N(i), 'uniform') ;
end

FigLogLog = figure(1);
loglog(h,err,'o-');
grid on
title({'-(cu'')''+au=f';'Dirichlet - Neumann'}); 
xlabel('err');
ylabel('hmax');

saveas(FigLogLog,'femDN.svg') ;

logh = log(h) ;
logerr = log(err) ; 
f=fit(logh,logerr,'poly1') ; 
FitLogLog = figure(2);
title({'-(cu'')''+au=f - loglogfit';'DN'})
plot(f,logh,logerr,'o') ;
title({'-(cu'')''+au=f';'Dirichlet - Neumann'}) ;  
xlabel('log(err)') ;
ylabel('log(hmax)') ;
saveas(FitLogLog, 'femDN-error.svg') ;

first=N(1)
last=N(lN)
disp(['boh ',num2str(first), ' Nmax ', num2str(last)])