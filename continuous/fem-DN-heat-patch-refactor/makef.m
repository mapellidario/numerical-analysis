function [fh] = makef(N, h, m, alpha, beta, gamma, BCtype)

switch strcat(BCtype)
    case 'DN'
        fh = zeros(N,1);
        for i=1:N-1
            fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
        end
        % first element for dirichlet in x=0
        fh(1) = fh(1) + alpha/h(1)*c(m(1)) - alpha*h(1)/4*rho(m(1));
        % last element for neumann in x=1
        fh(N) = h(N)/2*f(m(N)) + gamma;
    case 'DD'
        display('empty')
end


end