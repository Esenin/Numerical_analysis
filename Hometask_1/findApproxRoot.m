function [ x0 ] = findApproxRoot( poly, left, right, step)
    goodEps = 10^(-2);
    xL = left;    
    x0 = right + step;
    while (xL <= right - step)
        xR = xL + step;
        yLeft  = poly(xL);
        yRight = poly(xR);
        
        if (sign(yLeft) * sign(yRight) < 0 || abs((yLeft + yRight) / 2) < goodEps)
            x0 = (xL + xR) / 2;
            return;
        end;
        xL = xR;
    end;
end

