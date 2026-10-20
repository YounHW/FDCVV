function RUN_ONE_CASE(method, benchmark, repeat)
% Reproduce one complete 30-update trajectory in a fresh MATLAB session.
% This entry does not resume or overwrite an existing trajectory.
root=fileparts(mfilename('fullpath'));
config=readtable(fullfile(root,'configurations.csv'),'TextType','string');
q=config.method==string(method) & config.benchmark==string(benchmark);
assert(sum(q)==1,'Unknown method/benchmark combination'); c=config(q,:);
assert(isscalar(repeat)&&repeat==fix(repeat)&&repeat>=c.repeat_min&&repeat<=c.repeat_max);
originalPath=path;originalFolder=pwd;
cleanup=onCleanup(@()restore_session(originalPath,originalFolder));
SETUP_PACKAGE(method);
global ABL_B ABL_REP;ABL_B=find(strcmp({'GL1D','Peaks2D','H3_historical'},benchmark));ABL_REP=repeat;
dest=fullfile(root,'results',char(method),char(benchmark),sprintf('rep_%02d',repeat));
assert(~isfolder(dest),'Result folder already exists; preserve it and use a fresh package copy');
mkdir(dest);
rng(c.initial_seed,'twister');
[y,lb,ub,x,M]=call_benchmarkFunctions(c.benchmark_id);
x=scale_vector_to_unity(lb,ub,x);A=[lb';ub'];
assert(size(x,1)==c.initial_n);
seed=c.run_seed_base+repeat;
cd(dest);rng(seed,'twister');rng_before=rng;
save('launch_inputs.mat','x','y','lb','ub','A','seed','repeat','method','benchmark','rng_before');
if ismember(method,{'FDCVV','LC2015'})
    [models,final_errors,single_errors,time_data]=initial_metamodel( ...
        M,A,method,x,y,30,'MS',repeat,@Matern32_matrix);
    curve=single_errors{2}.error_data;
    save('OutputData.mat','models','final_errors','single_errors','time_data','-v7.3');
else
    output=chooseSamplingMethod(M,x,y,A,30,1,'MS',c.method_id);
    result=output{c.method_id};models=result{2}{1};
    curve=result{4}{1}{2}.error_data;
    save('OutputData.mat','output','-v7.3');
end
assert(numel(models)==30 && height(curve)==30 && all(isfinite(curve.Variables),'all'));
for t=1:30
    model=models{t};assert(model.m==size(x,1)+t);
    assert(all(isfinite(model.theta)) && all(isfinite(model.R),'all'));
    assert(all(model.X>=0 & model.X<=1,'all') && size(unique(model.X,'rows'),1)==model.m);
    yy=zeros(model.m,1);
    for j=1:model.m,yy(j)=M(lb'+model.X(j,:).*(ub-lb)');end
    assert(max(abs(yy-model.Y(:)))<1e-8);
    if ~strcmp(method,'LHS_rebuild')
        assert(isequal(model.X(1:size(x,1),:),x));
        if t>1,assert(isequal(model.X(1:end-1,:),models{t-1}.X));end
    end
end
writetable(curve,'native_convergence.csv');
test=load(fullfile(root,'tests',[char(benchmark) '_test.mat']));common=zeros(30,5);
for t=1:30
 pred=zeros(size(test.U,1),1);for j=1:numel(pred),pred(j)=models{t}.predict(test.U(j,:));end
 rmse=sqrt(mean((pred-test.truth).^2));span=max(pred)-min(pred);assert(span>0);
 common(t,:)=[t models{t}.m rmse rmse/span 1-sum((pred-test.truth).^2)/sum((test.truth-mean(test.truth)).^2)];
end
writetable(array2table(common,'VariableNames',{'step','samples','common_rmse','common_nrmse','common_r2'}),'common_convergence.csv');
accepted=struct('status','PASS','method',method,'benchmark',benchmark, ...
    'repeat',repeat,'seed',seed,'steps',30,'matlab',version);
save('accepted.mat','accepted');
fprintf('PASS %s %s repeat %d\n',method,benchmark,repeat);
end
function restore_session(p,d)
cd(d);path(p);
end
