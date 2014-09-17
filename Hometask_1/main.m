%imports
addpath('../support/');


polyExample = [9 0 -9 -36 1];

%    x3 - 4x + 2 = 0
poly11 = [1 0 -4 2];

%    3x6 + 172x5 + 125x4 -600x3 +125x2 + 172x + 3
poly6 = [3 172 125 -600 125 172 3];


tic
%NewtonRoots(poly11)
NewtonRoots(poly6)

toc


