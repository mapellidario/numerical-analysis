function [fig] = drawsol (uh, uh0, x, kmax, BCtype)
fig = figure(1);
switch strcat(BCtype)
    case 'NN'
        %
        % initial condition
        plot ([x], [uh0], 'ro-' )
        %
        % plot the solution in every time step
        hold on
        for k=1:kmax
            plot ([x], [uh(:,k)], 'b' )
        end
    case 'ND'
        disp('empty')
end