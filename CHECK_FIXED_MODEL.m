function CHECK_FIXED_MODEL(method)
% Isolated smoke: no hyperparameter optimization and no new response sampling.
assert(strcmp(method,'LC2015'));[root,upstream]=SETUP_PACKAGE(method);
dest=fullfile(root,'validation','portable');if ~isfolder(dest),mkdir(dest);end
cd(dest);rng(57022,'twister');
X=[0;.14;.31;.53;.77;1];Y=[0;.9;-.3;.6;-.7;.2];
model=OK_model(@Matern32_matrix,X,Y,'MS',3);
[mu,var]=model.predict(X);assert(max(abs(mu(:)-Y))<1e-10);
assert(all(isfinite(var)) && all(var>-1e-10));
C=linspace(.001,.999,999)';
if strcmp(method,'LC2015')
    [x,info]=lc2015_select(model,C);
    assert(info.score(info.index)==max(info.score));
    cv=lc2015_loo(model.R,model.Y);max_loo_error=0;
    for i=1:numel(Y)
        mask=true(numel(Y),1);mask(i)=false;
        fold=OK_model(@Matern32_matrix,X(mask),Y(mask),'MS',3);
        pred=fold.predict(X(i));
        max_loo_error=max(max_loo_error,abs((Y(i)-pred)-cv.error(i)));
    end
    assert(max_loo_error<1e-10);
    [xs,~,search]=lc_search(model,[0;1]);
    assert(search.converged && search.batches>=4 && search.batches<=16);
    assert(min(abs(X-xs))>1e-10);
    files={'OK_model','Matern32_matrix','lc2015_select','lc2015_loo','lc_original_predict','lc_search','LC2015_function','lhs_scaled'};
else,error('Use FDCVV or LC2015');end
assert(isfinite(x) && x>=0 && x<=1 && min(abs(X-x))>1e-10);
resolved=struct;
for i=1:numel(files)
    p=which(files{i});assert(startsWith(lower(p),lower(root))||startsWith(lower(p),lower(upstream)));resolved.(files{i})=p;
end
pathparts=strsplit(path,pathsep);
for i=1:numel(pathparts)
    p=pathparts{i};assert(isempty(p)||startsWith(lower(p),lower(root))||startsWith(lower(p),lower(matlabroot))||startsWith(lower(p),lower(upstream)), ...
        'Unexpected external MATLAB path');
end
report=struct('status','PASS','method',method,'matlab',version,'toolboxes',{ver}, ...
    'fixed_theta',3,'new_fits',0,'new_truth_evaluations',0,'selected_x',x,'resolved',resolved);
f=fopen([method '_smoke.json'],'w');fprintf(f,'%s',jsonencode(report));fclose(f);
fprintf('PASS portable %s fixed-model smoke; no fit or new truth evaluation\n',method);
end
