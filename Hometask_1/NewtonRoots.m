% @arg p -- coeffs of polynom
% @result -- roots vector
function [ roots ] = NewtonRoots(poly, aprxPoly, derivate)
    roots = [];
%const
    EPS = 10^(-5);    

% find root bounds
    [ downN, upperN, downP, upperP ]  = rootBounds(aprxPoly);
    format shortG;
    fprintf('bounds: [%i : %i]U[%i : %i]\n', downN, upperN, downP, upperP);
    positiveMinStep = (abs(upperP - downP))^(0.5) * 10^(-4);
    negativeMinStep = (abs(upperN - downN))^(0.5) * 10^(-4);
    fprintf('negative step = %i, positive step = %i\n', negativeMinStep, positiveMinStep);
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isNegativeInterval = true;
    curStep = negativeMinStep;
    curRightBarrier = upperN;
    x0 = downN;
    while (x0 < upperP && length(roots) < length(aprxPoly) - 1)
        % find approximation
        x0 = findApproxRoot(poly, x0, curRightBarrier, curStep);
        if (x0 <= curRightBarrier)
            fprintf('----------------------------------\n');
            fprintf('x0 = %d\n', x0);
            newRoot = runNewtonAtPoint(poly, derivate, x0);
                       
            if (length(roots) == 0)
                roots = [roots newRoot];
            elseif (abs(newRoot - roots(length(roots))) > 2 * EPS)
                roots = [roots newRoot];
            else 
                fprintf('This root is not actual;\n');
            end;
        end;
        x0 = x0 + curStep;
        if (x0 > upperN && isNegativeInterval)
            x0 = downP;
            curStep = positiveMinStep;
            curRightBarrier = upperP;
            isNegativeInterval = false;
        end;
        fprintf('\n');
    end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

