% modi normali di N pesi e N+1 molle (o N molle) appesi
%
% ascisse = autovalori
%
clear all
close all
%
% numero di pesi 
%
N = 5;
%
% estremo libero oppure fissato
%
libero = false;
% libero = true;
%
A = zeros(N+1,N);
A1 = A;
A2 = A;
%
A1(1:N,1:N) = +eye(N,N);
A2(2:N+1,1:N) = -eye(N,N);
A = A1+A2;
%
C = eye(N+1);
%
M = eye(N,N);
%
f = 0.01*ones(N,1);
%
if libero
    %
    % togliamo una molla se l'estremo e' libero
    %
    A = A(1:N,:);
    C = eye(N,N);
    %
end
%
% A and C matrix defined 

%
% modificare C(i,i) per cambiare le caratteristiche della i-esima molla
%
% C(2,2) = 100;
C(1,1) = 100 ;
%
K = A'*C*A;
%
[V,D] = eig(K,M); % Kx=lambda Mx
% eig returns two matrixes, 
%
% frequenze
% diag(D) is the vector containing the eigenvalues
%
omega = sqrt(diag(D));
%
% posizione di equilibrio
%
ueq = K\f;
%
xb = linspace(0,1,N+2);
xb = xb(2:N+1)';
%
Amp = 0.1;
% phase = 0
dt = 0.1;
Tmax = 20;
%
%j=0;
figure(1)
for t=0:dt:Tmax
    for i=1:N
        u = ueq + Amp*sin(omega(i)*t).*V(:,i);
        e = A*u;
        if libero
            plot([omega(i); omega(i)*ones(N,1)],1-[0; xb+u],'-O','LineWidth',2,'MarkerSize',10)
        else
            plot([omega(i); omega(i)*ones(N,1); omega(i)],1-[0; xb+u; 1],'-O','LineWidth',2,'MarkerSize',10)
        end
        if N>1
            axis([omega(1)-(omega(N)-omega(1))/N omega(N)+(omega(N)-omega(1))/N -0.1 1.1]);
        else
            axis([omega(1)-0.1 omega(1)+0.1 -0.1 1.1]);
        end
        grid on
        hold on
        plot(omega(i),1,'x','MarkerSize',15,'LineWidth',2)
        if ~libero
            plot(omega(i),0,'x','MarkerSize',15,'LineWidth',2)
        end
        title(sprintf('tempo %3.1f/%3.1f',t,Tmax))
        xlabel('frequenze proprie')
    end
    pause(0.01)
    hold off
end
