function cv = lc2015_loo(R,y)
% Single-fidelity ordinary Kriging, fixed correlation parameters.
% LC2015 supplement: fold GLS mean and REML scale, n-p-1 (p=1).
% R is the correlation matrix including any explicitly chosen nugget.
y=y(:); n=numel(y);
assert(n>2 && isequal(size(R),[n n]) && all(isfinite(y)));
assert(all(isfinite(R(:))) && norm(R-R','fro')<=1e-12*max(1,norm(R,'fro')));
L=chol(R,'lower'); P=L'\(L\eye(n)); s=P*ones(n,1); a=sum(s);
q=diag(P)-s.^2/a; assert(all(q>0),'LC: nonpositive precision diagonal');
centered=y-y(1); betaOffset=sum(P*centered)/a; residual=centered-betaOffset; b=P*residual;
rss=residual'*b; foldRSS=rss-b.^2./q;
tol=128*eps*max([abs(rss);abs(b.^2./q);realmin])*max(1,norm(P,inf)*norm(R,inf));
assert(all(foldRSS>=-tol),'LC: negative fold residual scale');
scale=max(foldRSS,0)/(n-2);
cv=struct('error',b./q,'variance',scale./q,'fold_scale',scale, ...
 'mean',y-b./q,'precision_diagonal',q,'full_reml_scale',rss/(n-1));
end

