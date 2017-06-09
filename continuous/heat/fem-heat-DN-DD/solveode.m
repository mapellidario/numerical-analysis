function [uh, uh0] = solveode (N, x, dt, kmax, Kh, Mh, fh, BCtype)

switch strcat(BCtype)
    case 'DN'
        % we now use implicit euler. we save every step
        uh = zeros(N,kmax) ;
        % uh(;,k) is the solution at step k.
        %
        % OSS: we do not have the initial condition in uh(:,0) !
        % we have to define the initial condition in a separate column vector
        
        uh0 = zeros(N,1) ;
        for i=1:N
            % interpolating initial condition u0
            uh0(i) = u0(x(i)) ;
            % uh0(i) = rand() ;
        end
    case 'DD'
        % we now use implicit euler. we save every step
        uh = zeros(N-1,kmax) ;
        % uh(;,k) is the solution at step k.
        %
        % OSS: we do not have the initial condition in uh(:,0) !
        % we have to define the initial condition in a separate column vector
        
        uh0 = zeros(N-1,1) ;
        for i=1:N-1
            % interpolating initial condition u0
            uh0(i) = u0(x(i)) ;
            % uh0(i) = rand() ;
        end
end


%
% horrible, but we have to pay that vectors start from 1 in matlab
% and we want to adhere to conventions used during lectures
% implicit euler step
for k=0:kmax-1
    % uh0 = uh(:,0) -> uh(:,1)
    if k==0
        uh(:,k+1) = (1/dt*Mh + Kh)\(1/dt*Mh*uh0+fh) ;
    else
        % uh0 = uh(:,k) -> uh(:,k+1)
        % step of implicit Euler
        uh(:,k+1) = (1/dt*Mh + Kh)\(1/dt*Mh*uh(:,k)+fh) ;
    end
end
%


end