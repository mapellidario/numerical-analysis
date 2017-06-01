function [] = drawpatch (N, kmax, x, uh, uh0, dt, alpha);

figure(2)
% draw the solution as a surface
% we use the patch command, which is supergeneric!
% being the grid cartesian we can use a more specific command, but
% this command is pretty generic and is nice to use it

% we assume that there is Neumann in x=1,
for i=1:N-1
    for k=1:kmax-1
        % point 1
        t1 = k*dt ;
        x1 = x(i) ;
        u1 = uh(i,k) ;
        % point 2
        t2 = (k+1)*dt ;
        x2 = x1 ;
        u2 = uh(i,k+1) ;
        % point3
        t3 = t2 ;
        x3 = x(i+1) ;
        u3 = uh(i+1,k+1) ;
        % point 4
        t4 = t1 ;
        x4 = x3 ;
        u4 = uh(i+1,k) ;
        %
        % uh in scala di colore
        % patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4]) ;
        % la voglio tutta bianca, quindi le varie u le considera come
        % altezza
        patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;

    end
end

% t=0
for i=1:N-1
    % point 1
    t1 = 0 ;
    x1 = x(i) ;
    u1 = uh0(i) ;
    % point 2
    t2 = dt ;
    x2 = x1 ;
    u2 = uh(i,1) ;
    % point3
    t3 = t2 ;
    x3 = x(i+1) ;
    u3 = uh(i+1,1) ;
    % point 4
    t4 = t1 ;
    x4 = x3 ;
    u4 = uh0(i+1) ;
    %
    % uh in scala di colore
    % patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4]) ;
    % la voglio tutta bianca, quindi le varie u le considera come
    % altezza
    patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;
end


% x=0
for k=1:kmax-1
    % point 1
    t1 = k*dt ;
    x1 = 0 ;
    u1 = alpha ;
    % point 2
    t2 = (k+1)*dt ;
    x2 = x1 ;
    u2 = alpha ;
    % point3
    t3 = t2 ;
    x3 = x(1) ;
    u3 = uh(1,k+1) ;
    % point 4
    t4 = t1 ;
    x4 = x3 ;
    u4 = uh(1,k) ;
    %
    % uh in scala di colore
    % patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4]) ;
    % la voglio tutta bianca, quindi le varie u le considera come
    % altezza
    patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;

end

% x=0, t=0
% point 1
t1 = 0 ;
x1 = 0 ;
u1 = alpha ;
% point 2
t2 = dt ;
x2 = x1 ;
u2 = alpha ;
% point3
t3 = t2 ;
x3 = x(1) ;
u3 = uh(1,1) ;
% point 4
t4 = t1 ;
x4 = x3 ;
u4 = uh0(1) ;
patch([t1 t2 t3 t4],[x1 x2 x3 x4], [u1 u2 u3 u4], 'w') ;


end