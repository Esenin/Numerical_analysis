% Newton method for equation system
function [x, y] = runNewtonAtPoint(f, g, fdx, fdy, gdx, gdy, x1, y1, a, k)
%   equation system:  { f(x,y) = 0.
%                     { g(x,y) = 0
% fdx, fdy \  
%           > derivatives
% gdx, gdy /

% x1, y1 -- start point
% x, y  -- arrays of iterative solutions
%   const:
EPS = 10^(-5);

    x = zeros(1, n);
    y = zeros(1, n);
    x(1) = x1;
    y(1) = y1;
    k = 1;
    while(k == 1 || abs(x(k) - x(k - 1)) > EPS || abs(y(k) - y(k - 1)) > EPS)
        xk = x(k);
        yk = y(k);
        a1 = f(xk, yk, a); b1 = fdx(xk, yk); c1 = fdy(xk, yk);
        a2 = g(xk, yk, k); b2 = gdx(xk, yk); c2 = gdy(xk, yk);
        det = b1 * c2 - c1 * b2;
        if (det ~= 0)
            detx = a1 * c2 - c1 * a2;
            dety = b1 * a2 - a1 * b2;
            x(k + 1) = xk - detx / det;
            y(k + 1) = yk - dety / det;
        else
            warning('Failed to execute the function Newton. Maybe incorrect x1, y1.');
            return;
        end;
        k = k + 1;
    end;
end