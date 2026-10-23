function values=COMMON_TEST_METRICS(pred,truth,truth_span)
% Columns: RMSE, legacy prediction-range NRMSE, R2, true range, primary NRMSE.
% RUN_ONE_CASE computes truth_span once from its fixed common test set.
pred=pred(:);truth=truth(:);
assert(numel(pred)==numel(truth) && numel(truth)>1);
assert(all(isfinite(pred)) && all(isfinite(truth)));
assert(isscalar(truth_span) && isfinite(truth_span) && truth_span>0);
sst=sum((truth-mean(truth)).^2);assert(sst>0);
sse=sum((pred-truth).^2);rmse=sqrt(sse/numel(truth));
prediction_span=max(pred)-min(pred);
legacy=NaN;
if prediction_span>0,legacy=rmse/prediction_span;end
values=[rmse legacy 1-sse/sst truth_span rmse/truth_span];
end
