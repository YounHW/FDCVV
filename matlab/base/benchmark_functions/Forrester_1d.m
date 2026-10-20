function [ y, lb, ub, x, M ] = Forrester_1d()
 


n = 1;
lb = 0;        % lower bound
ub = 1;        % upper bound

% Initial samples
x = scaled_TPLHD(5, lb, ub);  

M = @(xx) Forrester_function(xx);

y = zeros(size(x,1),1);
for i = 1:size(x,1)
    y(i,1) = M(x(i,:));
end

end


function [y] = Forrester_function(xx)

x1 = xx(1);

y = (6*x1 - 2)^2 * sin(12*x1 - 4);

end