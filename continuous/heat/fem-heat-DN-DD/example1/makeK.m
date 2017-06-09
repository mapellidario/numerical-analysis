function [Kh]=makeK(N, h, m, BCtype)

switch strcat(BCtype)
    case 'DN'
        % Kh = zeros(N,N);
        Kh = sparse(N,N);
        %
        for i=1:N-1
            %
            if i>1
                % only from second row
                Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
            end
            %
            Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
            %
            Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
        end
        % last row
        % Neumann homogenous
        Kh(N,N-1) = -1/h(N)*c(m(N)) ;
        Kh(N,N) = +1/h(N)*c(m(N)) ;
    case 'DD'
        Kh = sparse(N-1,N-1);
        %
        for i=1:N-1
            %
            if i>1
                % this term is not present in the first row
                Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
            end
            %
            Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
            %
            if i<N-1
                % this term is not present in the last row. meditate on this!
                Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
            end
        end
end

end