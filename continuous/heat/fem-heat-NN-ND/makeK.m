function [Kh]=makeK(N, h, m, BCtype)

switch strcat(BCtype)
    case 'NN'
        % Kh = zeros(N,N);
        Kh = sparse(N+1,N+1);
        %
        for i=2:N
            Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
            %
            Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
            %
            Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
        end
        % first row: Neumann
        Kh(1,1) = + 1/h(2)*c(m(2)) ;
        Kh(1,2) = - 1/h(2)*c(m(2)) ;
        % last row: Neumann homogenous
        Kh(N+1,N) = -1/h(N+1)*c(m(N+1)) ;
        Kh(N+1,N+1) = +1/h(N+1)*c(m(N+1)) ;
    case 'ND'
        Kh = sparse(N,N);
        %
        for i=2:N
          %
          Kh(i,i-1) = - 1/h(i)*c(m(i)) ;
          %  
          Kh(i,i) = + 1/h(i)*c(m(i)) + 1/h(i+1)*c(m(i+1)) ;
          %
          if i<N
            Kh(i,i+1) = - 1/h(i+1)*c(m(i+1)) ;
          end
        end
        % first row: Neumann
        Kh(1,1) = + 1/h(2)*c(m(2)) ;
        Kh(1,2) = - 1/h(2)*c(m(2)) ;

end

end
