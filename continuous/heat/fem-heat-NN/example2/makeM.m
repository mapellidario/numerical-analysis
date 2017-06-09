function [Mh] = makeM(N, h, m, BCtype)

switch strcat(BCtype)
    case 'NN'
        % Mh = zeros(N,N);
        Mh = sparse(N+1,N+1);
        for i=2:N
            Mh(i,i-1) = + h(i)/4*rho(m(i)) ;
            %
            Mh(i,i) = + h(i)/4*rho(m(i)) + h(i+1)/4*rho(m(i+1)) ;
            %
            Mh(i,i+1) = h(i+1)/4*rho(m(i+1)) ;
        end
        % first row
        % Neumann
        Mh(1,1)   = h(2)/4*rho(m(2)) ;
        Mh(1,2)   = h(2)/4*rho(m(2)) ;
        % last row
        % Neumann
        Mh(N+1,N)     = h(N+1)/4*rho(m(N+1)) ;
        Mh(N+1,N+1)   = h(N+1)/4*rho(m(N+1)) ;
        %
        % >> whos
        % >> spy(Mh)
    case 'ND'
        disp('empty')
end


end
