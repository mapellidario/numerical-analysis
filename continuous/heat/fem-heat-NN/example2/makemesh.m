function [x, h, m, N]=makemesh(meshtype, N)

    switch strcat(meshtype)
        case 'uniform'
            % we start with a uniform mesh (x_uniform)
            % >> linspace(xmin, xmax, M)
            % divide the interval [xmin,xmax] in M-1 intervals (with M points)
            %
            xu = linspace(0,1,N+1);
            %
            % problem: x_1=0, x_N+1 = 1
            %
            % xu = xu(2:end);
            %
            x = xu ;
            %
            % this new x is the old x, but starts from 2nd element
            % this aligns matlab vector with the one used during lectures
            % in matlab vectors strart from 1, there is no index = 0 in arrays.
        case 'transformed'
            % contruct non uniform mesh, by transforming the current mesh
            % using the point operator for term by term operations
            % >> x = [1 2 3 4]
            % >> x.^2;
            % x: [1 4 9 16]
            % x = xu.^2 ;
            %
            disp('not implemented yet')
            %
        case 'random'
            % random mesh, N random points in (0,1)
            %
            x = rand(1,N+1) ;
            % sort the array, from lower to higher
            x = sort(x) ;
            % add x=0 and x=1
            x = [0 x 1] ;
            % remove repeated elements
            x = unique(x) ;
            %
            % remove x=0
            x = x(2:end) ;
            % redefine N
            N = length(x) ;
        otherwise
            %
            disp('wrong mesh name')
            %
    end
    %
    % compute h and m
    % this works for a generic mesh
    %
    h = zeros(1,N+1);
    m = zeros(1,N+1);

    %
    % i-1 = 0 for i=1, so this cycle starts from 2
    %
    % changed!
    % h(1) = x(1) - 0;
    % m(1) = (x(1) + 0)/2;
    % ACHTUNG h(1) and m(1) are not used anywhere!
    % just initialise them
    h(1)= 0. ;
    m(1)= 0 ;
    
    for i=2:N+1
        h(i) = x(i)-x(i-1);
        m(i) = (x(i-1)+x(i))/2 ;
    end
    %
    % inspect the mesh
    % >> bar(h)
    %

end
