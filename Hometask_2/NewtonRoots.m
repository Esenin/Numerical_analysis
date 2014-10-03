% @result -- roots vector
function [ roots ] = NewtonRoots(f, g, fdx, fdy, gdx, gdy, a, k)
    roots = [];
%const
    EPS = 10^(-5);
    APRXEPS = 10^(-2);    

% find root bounds
    startX = 0;
    stopX = 2;
    startY = 0;
    stopY = 2;
    
    format shortG;
    fprintf('bounds: X=[%i : %i], Y=[%i : %i]\n', startX, stopX, startY, stopY);
    aprxStep = .01;
    
    fprintf('aprxStep = %i\n', aprxStep);
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x0 = startX;
    y0 = startY;
    % find approximation
    for ai = a
    	for ki = k
            printf('------------------------------------------');
            printf('----    a = %i,     k = %i \n', ai, ki);
            while (x0 <= stopX)
                if (x^2 > k)
                    x0 = x0 + aprxStep;
                    continue;
                end;
                y01 = (k - x^2)^(0.5);
                y02 = -1 * (k - x^2)^(0.5);
                
                if (abs(f(x0, y01, ai)) < APRXEPS && abs(g(x0, y01, ki)) < APRXEPS)
                    % near root point
                    [xks, yks] = runNewtonAtPoint(f, g, fdx, fdy, gdx, gdy, x0, y01, ai, ki);
                end;
                if (abs(f(x0, y02, ai)) < APRXEPS && abs(g(x0, y02, ki)) < APRXEPS)
                    % near root point with negative Y 
                    [xks, yks] = runNewtonAtPoint(f, g, fdx, fdy, gdx, gdy, x0, y02, ai, ki);
                end;
                
            
                x0 = x0 + aprxStep;
            end;
            
            
        end;
    end;
                        
                        
    
    
 
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

