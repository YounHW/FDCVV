function [xnew,info] = lc2015_select(model,C,duplicateTolerance)
% LC2015 Eq.(8)-(9), s=1; R reference voronoi()/adjustment().
% Finite supplied candidate pool; first-index deterministic tie breaking.
% Adapter accepts RevisionKriging with X,Y,R,predict.
if nargin<3, duplicateTolerance=1e-10; end
assert(size(C,2)==size(model.X,2) && all(isfinite(C(:))));
cv=lc2015_loo(model.R,model.Y);
% Use latent-function held-out variance, consistent with model.predict.
jitter=0;
if isprop(model,'diagnostics'), jitter=model.diagnostics.jitter; end
cv.variance=cv.variance-jitter*cv.fold_scale;
assert(all(cv.variance>=0),'LC: negative latent fold variance');
ratio=zeros(size(cv.error)); positive=cv.variance>0;
ratio(positive)=cv.error(positive).^2./cv.variance(positive);
assert(~any(~positive & cv.error~=0),'LC: undefined nonzero-error/zero-variance ratio');
bestDistance=inf(size(C,1),1); owner=ones(size(C,1),1);
for i=1:size(model.X,1)
 ds=sum((C-model.X(i,:)).^2,2); take=ds<bestDistance;
 bestDistance(take)=ds(take); owner(take)=i;
end
[~,variance]=lc_original_predict(model,C); score=variance(:).*(1+ratio(owner));
valid=bestDistance>duplicateTolerance^2; score(~valid)=-inf;
assert(any(valid) && all(isfinite(score(valid))),'LC: no finite valid candidate');
if all(score(valid)==0)
 % Explicit constant-response convention: maximin exploration.
 fallback=bestDistance; fallback(~valid)=-inf; [~,at]=max(fallback);
else
 [~,at]=max(score);
end
xnew=C(at,:);
info=struct('loo',cv,'ratio',ratio,'owner',owner,'score',score, ...
 'index',at,'jitter',jitter,'search','supplied candidate pool');
end
