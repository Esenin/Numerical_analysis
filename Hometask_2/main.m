clear 
clc

%imports
addpath('../support/');
format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  { tg(ax + y) - axy = 0.3
%  { x^2 + y^2 = k
%    a = -0.6 : -0.1 : -1.1
%    k = 1.3 : 0.1 : 1.8
a = -0.6 : -0.1 : -1.1;
k = 1.3 : 0.1 : 1.8

f = @(x, y, a) tan(a*x + y) - a*x*y - 0.3;
g = @(x, y, k) x^2 + y^2 - k;

fdx = @(x, ~) 2 * x; 
fdy = @(~, y) 2 * y;
gdx = @(x, y) 3 * x ^ 2 * y; 
gdy = @(x, ~) x ^ 3; 



roots = NewtonRoots(f, g, fdx, fdy, gdx, gdy, a, k);


%af = 0 : 0.01 : 2 * pi; 
%rf = sqrt(2) * ones(1, size(af, 2)); 
%xg = -3 : 0.01 : 3;
%yg = xg .^ (-3);

hold on;
grid on;
n = length(x);
graph1 = polar(af, rf, 'k');
set(graph1, 'LineWidth', 2);
graph2 = 


%plot (xg, yg, 'k--', ...
%      x(n), y(n), 'bp', ...
%      x1, y1, 'ro', ...
%      'LineWidth', 2);
%  
axis([0, 3, 0, 3]);

legend('x^2 + y^2 = 2', ...
       'x^3 y = 1', ...
       'unaccurate solution', ...
       'sample text!!11');

tic
toc


