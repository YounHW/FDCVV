function [ y, lb, ub, x, M ] = GoldsteinPrice_2d()
 


n = 2;
lb = [-2; -2];        % lower bound
ub = [ 2;  2];        % upper bound

% Initial samples
x = scaled_TPLHD(5, lb, ub);  

M = @(xx) GoldsteinPrice_function(xx);

y = zeros(size(x,1),1);
for i = 1:size(x,1)
    y(i,1) = M(x(i,:));
end

end


function [y] = GoldsteinPrice_function(xx)

x1 = xx(1);
x2 = xx(2);

term1 = 1 + (x1 + x2 + 1)^2 * ...
    (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2);

term2 = 30 + (2*x1 - 3*x2)^2 * ...
    (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2);

y = term1 * term2;

end