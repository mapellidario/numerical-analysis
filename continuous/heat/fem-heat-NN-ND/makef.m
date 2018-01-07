function [fh] = makef(N, h, m, BCtype, delta, gamma, beta)

% notice that the fh vector has different sizes acording to the BC!
% this is why, even though for i=1:N-1 the vector is the same, it is
% inside the switch nonetheless.

switch strcat(BCtype)
    case 'NN'
        fh = zeros(N+1,1);
        for i=2:N
            fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
        end
        % first element for neumann non homogenous in x=0
        % considered K and M matrix effect
        fh(1) = h(2)/2*f(m(2)) + delta ;
        % last element for neumann non homogenous in x=1
        fh(N+1) = h(N+1)/2*f(m(N+1)) + gamma ;
    case 'ND'
        fh = zeros(N,1);
        for i=2:N
            fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
        end
        % first element for neumann non homogenous in x=0
        % considered K and M matrix effect
        fh(1) = h(2)/2*f(m(2)) + delta ;
        % last element for dirichlet non homogenous in x=1
        fh(N) = fh(N) + beta/h(N+1)*c(m(N+1)) - beta*h(N+1)/4*rho(m(N+1)) ;
end


end
