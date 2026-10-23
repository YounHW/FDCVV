function CHECK_COMMON_METRICS
% Check the reporting metric without fitting models or evaluating new samples.
root=fileparts(mfilename('fullpath'));
truth=[0;2;4];span=4;
a=COMMON_TEST_METRICS([1;2;3],truth,span);
b=COMMON_TEST_METRICS([-1;2;5],truth,span);
assert(abs(a(1)-sqrt(2/3))<1e-14 && a(1)==b(1));
assert(a(5)==b(5) && abs(a(5)-a(1)/4)<1e-14);
assert(a(2)~=b(2));
constant=COMMON_TEST_METRICS([2;2;2],truth,span);
assert(isnan(constant(2)) && isfinite(constant(5)));
perfect=COMMON_TEST_METRICS(truth,truth,span);
assert(perfect(1)==0 && perfect(5)==0 && perfect(3)==1);
benchmarks={'GL1D','Peaks2D','H3_historical'};
for i=1:3
    test=load(fullfile(root,'tests',[benchmarks{i} '_test.mat']));
    truth=test.truth(:);span=max(truth)-min(truth);
    a=COMMON_TEST_METRICS(truth+.01*span,truth,span);
    assert(abs(a(5)-.01)<1e-12 && a(4)==span);
end
fprintf('PASS: fixed truth-range metric, legacy field, constant/perfect predictions and three common test sets.\n');
end
