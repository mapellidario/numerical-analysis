clear all
close all
%
% definiamo la suddivisione
%
% cominciamo con una uniforme
%
% linspace(a,b,M) divide 
% l'intervallo [a,b] in 
% M-1 intervallini
%
% noi vogliamo N intervallini
%
N = 2;
%
x = linspace(0,1,N+1);
%
x = x(2:end);
%
% x è quello delle note
% ma non c'è x(0) !!!!
%
% calcoliamo gli h e gli m
%
h = zeros(1,N);
m = zeros(1,N);
%
h(1) = x(1)-0;
m(1) = (0+x(1))/2;
%
for i=2:N
    h(i) = x(i)-x(i-1);
    m(i) = (x(i-1)+x(i))/2;
end
%
Kh = zeros(N,N);
%
for i=1:N-1
    %
    if i>1
        % prima riga!!!
        Kh(i,i-1) = -1/h(i)*c(m(i));
    end
    %
    Kh(i,i) = +1/h(i)*c(m(i)) + ...
              +1/h(i+1)*c(m(i+1));
    %
    Kh(i,i+1) = -1/h(i+1)*c(m(i+1));
    %
end
% ultima riga
Kh(N,N-1) = -1/h(N)*c(m(N));
Kh(N,N) = +1/h(N)*c(m(N));
%
% termine noto
%
fh = zeros(N,1);
%
for i=1:N-1
    fh(i) = h(i)/2*f(m(i)) + ...
            h(i+1)/2*f(m(i+1));
end
% ultimo elemento:
fh(N) = h(N)/2*f(m(N));
%
uh = Kh\fh; % uh e' un vettore colonna
%
plot([0 x],[0; uh],'k-')




























    
























