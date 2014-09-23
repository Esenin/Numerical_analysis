function [ root ] = runNewtonAtPoint( poly, derivate,  x0 )
%const
    EPS = 10^(-5);    
% make approx polynom
    counter = 0;
    xk = x0;
    xkn = x0 + EPS + 1;
    
    fprintf('\ti =\t|\tX =\t\t\t\t|\tY =\t\n');
    while (abs(xk - xkn) >= EPS)
        counter = counter + 1;
                
        xkn = xk - poly(xk) / derivate(xk);
        fprintf('\t%i\t|\t%i\t|\t%i\t\n', counter, xk, poly(xkn));
        temp = xk;
        xk = xkn;
        xkn = temp;
    end;
    
    root = xk;
end

