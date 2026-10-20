function [ y, lb, ub,x, M ] = Mod_Gramacy_Lee_1d()
 


n = 1;
lb = -1.5*ones(n,1);        % lower bound
ub = 6*ones(n,1); 


% Initial samples
x = scaled_TPLHD(5,lb,ub);  

M = @(xx) Gramacy_Lee_function(xx);

y = zeros(size(x,1),1);
for i=1:size(x,1)
    y(i,1) = M(x(i,:));
end



end


function y = Gramacy_Lee_function(xx)

% 只调整这些常数：频率越小，波越疏
f_left  = 3.5;   % 原来是 6
f_right = 2.5;   % 原来是 4

%term1 = 50*sin( f_left *pi*xx ) / (cos(xx));   % 原来 sin(6*pi*xx)
den = cos(xx);
delta = 1e-2;  % 可试 1e-3 ~ 1e-1
if abs(den) < delta
    den = sign(den)*delta + (den==0)*delta;
end
term1 = 50*sin( f_left*pi*xx) / den;
term2 = (xx-1)^4;

if xx>0.5 && xx<=2.5
   term2 = 0.8*(xx-1)^4;
   term1 = 20*sin( f_left *pi*xx)/(3*(xx));    % 原来 sin(6*pi*xx)
   epsilon = -100.0;
elseif xx>2.5
   term2 = 0.5*(xx-1)^3;
   term1 = 100*sin( f_right*pi*xx)/((xx));   % 原来 sin(4*pi*xx)
   epsilon = -50.0;
else
   epsilon= -30.0;
end
y = term1 + term2+epsilon;
end


