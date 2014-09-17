function [ x0 ] = findApproxRoot( polyn, left, right )
    [left right]

% step 
    step = abs(right - left) / 1000;
    goodEps = 0.01;
    x = left : step : right;
    iters = length(x);
    
    x0 = right + step;
    for i = 1 : iters
        %calcPoly(polyn, x(i))
        if (abs(calcPoly(polyn, x(i))) < goodEps)
            x0 = x(i);
            return;
        end;
    end;
end

