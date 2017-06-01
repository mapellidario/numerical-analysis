function [] = drawsol (uh, uh0, x, kmax, BCtype, alpha, beta)

switch strcat(BCtype)
    case 'DN'
        figure(1)
        %
        % initial condition
        plot ([0 x], [alpha;uh0], 'ro-' )
        %
        % plot the solution in every time step
        hold on
        for k=1:kmax
            plot ([0 x], [alpha; uh(:,k)], 'b' )
        end
    case 'DD'
        display('empty')
end