function [Mh] = makeM(N, h, m, BCtype)

switch strcat(BCtype)
    case 'DN'
        % Mh = zeros(N,N);
        Mh = sparse(N,N);
        for i=1:N-1
            %
            if i>1
                % only from second row
                Mh(i,i-1) = + h(i)/4*rho(m(i)) ;
            end
            %
            Mh(i,i) = + h(i)/4*rho(m(i)) + h(i+1)/4*rho(m(i+1)) ;
            %
            Mh(i,i+1) = h(i+1)/4*rho(m(i+1)) ;
        end
        % last row
        % Neumann homogenous
        Mh(N,N-1) = h(N)/4*rho(m(N)) ;
        Mh(N,N)   = h(N)/4*rho(m(N)) ;
        %
        % >> whos
        % >> spy(Mh)
    case 'DD'
        display('empty')
        % I still have to compute this!
end


end
