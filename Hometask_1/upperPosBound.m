% @arg c - coeff vector of polynom with zeros. First value - max degree,
% a_0 > 0
% @result upper bound of positive roots
function [ upP ] = upperPosBound( c )
    % find upper bounds
    firstNegCoefIndex = -1;
    maxNegate = 1;
    for i = 1 : length(c)
        if (firstNegCoefIndex < 0 && c(i) < 0)
            firstNegCoefIndex = i - 1;
        end;
        if (c(i) < maxNegate)
            maxNegate = c(i);
        end;
    end;
    maxNegate = abs(maxNegate);
    
    
   % disp(['subRootExpr:']);
   % maxNegate / c(1)
   % disp(['degree of root:']);
  %  (1 / firstNegCoefIndex)
    
    upP = 1 + power(maxNegate / abs(c(1)), (1 / firstNegCoefIndex));
    
end

