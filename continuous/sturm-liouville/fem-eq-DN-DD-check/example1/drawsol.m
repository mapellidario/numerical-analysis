function [fig] = drawsol (uh, x, BCtype, alpha, beta)

fig = figure(1) ;

switch strcat(BCtype)
    case 'DN'
        plot ([0 x], [alpha; uh], 'o-' )
    case 'DD'
        plot([0 x], [alpha; uh; beta], 'o-');
end