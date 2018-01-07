function [fig] = drawsol (uh, uh0, x, kmax, BCtype, alpha, beta)
fig = figure(1) ;
switch strcat(BCtype)
    case 'DN'
        %
        % initial condition
        plot ([0 x], [alpha;uh0], 'ro-' )
        %
        % plot the solution in every time step
        hold on
        for k=1:kmax
            % uh(:,k) is the kth column, append alpha at the beginning
            %plot (x,y,options)
            plot ([0 x], [alpha; uh(:,k)], 'b' )
        end
    case 'DD'
        %
        % initial condition
        plot([0 x], [alpha; uh0; beta], 'ro-');
        %
        % plot the solution in every time step
        hold on
        for k=1:kmax
            plot ([0 x], [alpha; uh(:,k); beta], 'b' )
        end
end