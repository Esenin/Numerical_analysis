% @arg p -- coeffs of polynom
% @result -- roots vector
function [ roots ] = NewtonRoots( p )
    roots = {};
%const
    eps = 10^(-5);

% find root bounds
    [ downN, upperN, downP, upperP ]  = rootBounds(p)
    positiveMinStep = abs(upperP - downP) / 100
    negativeMinStep = abs(upperN - downN) / 100
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isNegativeInterval = true;
    curStep = negativeMinStep;
    curRightBarrier = upperN;
    x0 = downN;
    while (x0 < upperP && length(roots) < length(p) - 1)
        
        % find approximation
        x0 = findApproxRoot(p, x0, curRightBarrier);
        % make approx polynom
        apprxP = @(x) calcPoly(p, x0) + (x - x0) * calcPoly(getDerivate(p), x0);
        apprxx = apprxP(x0);

        xk = x0;
        xkn = x0 + eps + 1;
        while (abs(xk - xkn) > eps)
            xkn = xk - (calcPoly(p, xk) / calcPoly(getDerivate(p), xk));
            temp = xk;
            xk = xkn;
            xkn = temp;
        end;
        roots = [roots xk];
        x0 = x0 + curStep;
        if (x0 > upperN && isNegativeInterval)
            x0 = downP;
            curStep = positiveMinStep;
            curRightBarrier = upperP;
            isNegativeInterval = false;
        end;
    end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

