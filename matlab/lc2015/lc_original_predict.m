function [mu,variance] = lc_original_predict(model,C)
% Algebraically identical to original OK_model.predict, batched for LC search.
% Original product Matern32 correlation, GLS mean and ML process variance.
assert(strcmp(func2str(model.auto_correlation_function),'Matern32_matrix'));
beta=model.compute_beta_hat(model.R);
scale=model.compute_sigma_sq_hat(model.R,beta);
weights=model.R\(model.Y-model.F*beta);
rf=model.R\model.F; denom=model.F'*rf;
mu=zeros(size(C,1),1); variance=mu;
for first=1:2048:size(C,1)
    ids=first:min(first+2047,size(C,1));
    r=ones(model.m,numel(ids));
    for d=1:size(C,2)
        z=sqrt(3)*model.theta(d)*abs(model.X(:,d)-C(ids,d)');
        r=r.*((1+z).*exp(-z));
    end
    solved=model.R\r; u=model.F'*solved-1;
    mu(ids)=beta+r'*weights;
    variance(ids)=scale*(1-sum(r.*solved,1)+u.^2/denom);
end
end
