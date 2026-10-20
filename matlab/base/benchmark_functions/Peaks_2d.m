function [ y, lb, ub, x, M ] = Peaks_2d()
 


n = 2;
lb = [-4; -4];        % lower bound
ub = [ 4;  4];        % upper bound

% Initial samples
x = scaled_TPLHD(10, lb, ub);  

M = @(xx) Peaks_function(xx);

y = zeros(size(x,1),1);
for i = 1:size(x,1)
    y(i,1) = M(x(i,:));
end

end


function [y] = Peaks_function(xx)

x1 = xx(1);
x2 = xx(2);

term1 = 3*(1 - x1)^2 * exp(-(x1^2) - (x2 + 1)^2);
term2 = -10*(x1/5 - x1^3 - x2^5) * exp(-x1^2 - x2^2);
term3 = -(1/3) * exp(-(x1 + 1)^2 - x2^2);

y = term1 + term2 + term3;

end