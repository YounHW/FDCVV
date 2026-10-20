function [ y, lb, ub, x, M ] = MultiPeak_1d()
 


n = 1;
lb = -4;
ub = 4;

% Initial samples
x = scaled_TPLHD(3, lb, ub);  

M = @(xx) MultiPeak_function(xx);

y = zeros(size(x,1),1);
for i = 1:size(x,1)
    y(i,1) = M(x(i,:));
end

end


function [y] = MultiPeak_function(xx)

x1 = xx(1);

y = 3*(1 - x1)^2 * exp(-x1^2) ...
    - 10*(x1/5 - x1^3) * exp(-x1^2) ...
    - 1/3 * exp(-(x1 + 1)^2);

end