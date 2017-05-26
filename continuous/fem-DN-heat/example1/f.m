function y = f(x)
%    y = 1 ;
%    y = 2*x*(5*sin(5*x) - 1/(x + 1)) + sin(x)^2*(cos(5*x) + log(x + 1)) + (25*cos(5*x) + 1/(x + 1)^2)*(x^2 + 1) ;
    y = sin(x)^2*(log(x + 1) + sin(5*x)) - 2*x*(5*cos(5*x) + 1/(x + 1)) + (x^2 + 1)*(25*sin(5*x) + 1/(x + 1)^2) ;

%    y = sin(6*x)*x^2 ;
%    y = (x^2 + 1)*(25*sin(5*x) + 1/(x + 1)^2) - 2*x*(5*cos(5*x) + 1/(x + 1));

end
