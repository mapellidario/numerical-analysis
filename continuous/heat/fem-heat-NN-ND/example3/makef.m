function [fh] = makef(N, h, m, BCtype, delta, gamma)

% notice that the fh vector has different sizes acording to the BC!
% this is why, even though for i=1:N-1 the vector is the same, it is 
% inside the switch nonetheless.

switch strcat(BCtype)
    case 'NN'
        fh = zeros(N+1,1);
        for i=2:N
            fh(i) = h(i)/2*f(m(i)) + h(i+1)/2*f(m(i+1)) ;
        end
        % first element for neumann non homogenous in x=0
        % considered K and M matrix effect
        fh(1) = h(2)/2*f(m(2)) + delta ;
        % last element for neumann non homogenous in x=1
        fh(N+1) = h(N+1)/2*f(m(N+1)) + gamma ;
    case 'ND'
        disp('empty')
end


end