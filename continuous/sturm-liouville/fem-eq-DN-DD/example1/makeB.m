function [Bh] = makeB (N, h, m, BCtype) ;

switch strcat(BCtype)
    case 'DN'
        % Mh = zeros(N,N);
        Bh = sparse(N,N);
        for i=1:N-1
            %
            if i>1
                % only from second row
                % Bh(i,i-1) = + b(m(i))/2 ;
                Bh(i,i-1) = - b(m(i))/2 ;
            end
            %
            Bh(i,i) = + b(m(i))/2 - b(m(i+1))/2 ;
            % Bh(i,i) = - b(m(i))/2 + b(m(i+1))/2 ;
            %
            % Bh(i,i+1) = - b(m(i+1))/2 ;
            Bh(i,i+1) = + b(m(i+1))/2 ;
        end
        % last row
        % Neumann homogenous
        Bh(N,N-1) = - b(m(N))/2 ;
        Bh(N,N)   = + b(m(N))/2 ;
        %
        % >> whos
        % >> spy(Mh)
    case 'DD'
        Bh = sparse(N-1,N-1);
        for i=1:N-1
            %
            if i>1
                % only from second row
                % Bh(i,i-1) = b(m(i))/2 ;
                Bh(i,i-1) = - b(m(i))/2 ;
            end
            %
            Bh(i,i) = + b(m(i))/2 - b(m(i+1))/2 ;
            % Bh(i,i) = - b(m(i))/2 + b(m(i+1))/2 ;
            %
            if i<N-1
                % Bh(i,i+1) = - b(m(i+1))/2 ;
                Bh(i,i+1) = + b(m(i+1))/2 ;
            end
        end
end


end