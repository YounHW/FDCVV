function [ y, lb, ub, x, M ] = Ackley_3d()
 


n = 3;
lb = -32.768 * ones(n,1);      % lower bound
ub =  32.768 * ones(n,1);      % upper bound

% Initial samples
x = scaled_TPLHD(10, lb, ub);  

M = @(xx) Ackley_function(xx);

y = zeros(size(x,1),1);
for i = 1:size(x,1)
    y(i,1) = M(x(i,:));
end

end


function [y] = Ackley_function(xx)

x1 = xx(1);
x2 = xx(2);
x3 = xx(3);

n = 3;

a = 20;
b = 0.2;
c = 2*pi;

sum1 = x1^2 + x2^2 + x3^2;
sum2 = cos(c*x1) + cos(c*x2) + cos(c*x3);

term1 = -a * exp(-b * sqrt(sum1 / n));
term2 = -exp(sum2 / n);

y = term1 + term2 + a + exp(1);

end