% @result -- roots vector
function [ ] = NewtonRoots(f, g, fdx, fdy, gdx, gdy, a, k)
xks = zeros(1, 3);
yks = zeros(1, 3);

%const
    EPS = 10^(-5);
    APRXEPS = 0.06;    

% find root bounds
    startX = -2;
    stopX = 2;
    
    format shortG;
    fprintf('bounds: X=[%i : %i], Y=[-1 * (k - x^2)^(0.5) : (k - x^2)^(0.5)]\n', startX, stopX);
    aprxStep = .01;
    
    fprintf('aprxStep = %i\n', aprxStep);
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % find approximation
    for ai = a
        for ki = k
            fprintf('---------------------------------------------------------------\n');
            fprintf('----    a = %i,     k = %i \n', ai, ki);
            x0 = startX;
            lastX = -1 * 10^5;
            hasRoot = false;
            fA = @(x,y) f(x,y,ai);
            gK = @(x,y) g(x,y,ki);
            fdxA = @(x,y) fdx(x, y, ai);
            fdyA = @(x,y) fdy(x,y,ai);
            gdxK = @(x,y) gdx(x,y,ki);
            gdyK = @(x,y) gdy(x,y,ki);

            while (x0 <= stopX)
                if (x0^2 > ki)
                    x0 = x0 + aprxStep;
                    continue;
                end;
                y01 = (ki - x0^2)^(0.5);
                y02 = -1 * (ki - x0^2)^(0.5);
                
                if (abs(f(x0, y01, ai)) < APRXEPS && abs(g(x0, y01, ki)) < APRXEPS)
                    % near root point
                    [xks, yks] = runNewtonAtPoint(fA, gK, fdxA, fdyA, gdxK, gdyK, x0, y01);
                    if (abs(xks(end) - lastX) > EPS)
                        writeoutSeq(xks, yks, fA, gK);
                        hasRoot = true;
                    end;
                    lastX = xks(end);
                end;
                if (abs(f(x0, y02, ai)) < APRXEPS && abs(g(x0, y02, ki)) < APRXEPS)
                    % near root point with negative Y 
                    [xks, yks] = runNewtonAtPoint(fA, gK, fdxA, fdyA, gdxK, gdyK, x0, y02);
                    if (abs(xks(end) - lastX) > EPS)
                        writeoutSeq(xks, yks, fA, gK);
                        hasRoot = true;
                    end;
                    lastX = xks(end);
                end;
                
                x0 = x0 + aprxStep;
            end;
            if (~hasRoot)
                fprintf('----    Roots were not found..\n');
            end;
            
            fprintf('\n\n');
            
        end;
    end;
         
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

