function [fh] = makef(N, h, m, BCtype, alpha, beta, gamma)

% notice that the fh vector has different sizes acording to the BC!
% this is why, even though for i=1:N-1 the vector is the same, it is 
% inside the switch nonetheless.

switch strcat(BCtype)
    case 'DN'
        fh = zeros(N,1);
        for i=1:N-1
            fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
        end
        % first element for dirichlet non homogenous in x=0
        % considered K and M matrix effect
        fh(1) = fh(1) + alpha/h(1)*c(m(1)) + alpha/2*b(m(1)) - alpha*h(1)/4*a(m(1)) ;
        % last element for neumann non homogenous in x=1
        fh(N) = h(N)/2*f(m(N)) + gamma ;
    case 'DD'
        % 
        fh = zeros(N-1,1);
        % this part is common
        for i=1:N-1
            fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
        end
        % first element for dirichlet non homogenous in x=0
        % considered K and M matrix effect
        fh(1)   = fh(1)   + alpha/h(1)*c(m(1)) + alpha/2*b(m(1)) - alpha*h(1)/4*a(m(1)) ;
        % last element for neumann non homogenous in x=1
        fh(N-1) = fh(N-1) + beta/h(N)*c(m(N))  - beta/2*b(m(N))  - beta*h(N)/4*a(m(N))  ;
end


end