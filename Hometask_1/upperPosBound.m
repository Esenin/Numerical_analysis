% @arg c - coeff vector of polynom with zeros. First value - max degree,
% a_0 > 0
% @result upper bound of positive roots
function [ upP ] = upperPosBound( aprxPoly )
    isFirstCoeffNeg = aprxPoly(1) < 0;
    if (isFirstCoeffNeg)
        aprxPoly = -aprxPoly;
    end;
    % find upper bounds
    firstNegCoefIndex = -1;
    maxNegate = 1;
    for i = 1 : length(aprxPoly)
        if (firstNegCoefIndex < 0 && aprxPoly(i) < 0)
            firstNegCoefIndex = i - 1;
        end;
        if (aprxPoly(i) < maxNegate)
            maxNegate = aprxPoly(i);
        end;
    end;
    maxNegate = abs(maxNegate);
        
    upP = 1 + power(maxNegate / aprxPoly(1), (1 / firstNegCoefIndex));
   
end

