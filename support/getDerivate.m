function [ p1 ] = getDerivate( p )
% p - coeff of polynom
    n = length(p);
    p1 = zeros(n - 1, 1);
    
    for i = 1 : n - 1
        degree = n - i;
        p1(i) = p(i) * degree;
    end;
end

