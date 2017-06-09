function [Ah] = makeA (N, h, m, BCtype)

switch strcat(BCtype)
    case 'DN'
        % Ah = zeros(N,N);
        Ah = sparse(N,N);
        for i=1:N-1
            %
            if i>1
                % only from second row
                Ah(i,i-1) = + h(i)/4* a (m(i)) ;
            end
            %
            Ah(i,i) = + h(i)/4*a(m(i)) + h(i+1)/4*a(m(i+1)) ;
            %
            Ah(i,i+1) = h(i+1)/4*a(m(i+1)) ;
        end
        % last row
        % Neumann homogenous
        Ah(N,N-1) = h(N)/4*a(m(N)) ;
        Ah(N,N)   = h(N)/4*a(m(N)) ;
        %
        % >> whos
        % >> spy(Ah)
    case 'DD'
        Ah = sparse(N-1,N-1);
        for i=1:N-1
            %
            if i>1
                % only from second row
                Ah(i,i-1) = + h(i)/4*a(m(i)) ;
            end
            %
            Ah(i,i) = + h(i)/4*a(m(i)) + h(i+1)/4*a(m(i+1)) ;
            %
            if i<N-1
                Ah(i,i+1) = h(i+1)/4*a(m(i+1)) ;
            end
        end
end


end
