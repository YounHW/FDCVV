function CHECK_FDCVV_REPLAY
root=SETUP_PACKAGE('FDCVV');global ABL_B ABL_REP;
dest=fullfile(root,'validation','replay');if ~isfolder(dest),mkdir(dest);end;cd(dest);
f=load(fullfile(root,'tests','replay_fixture.mat'));rows=[];
for i=1:numel(f.fixtures)
 z=f.fixtures{i};ABL_B=z.b;ABL_REP=z.rep;
 model=OK_model(@Matern32_matrix,z.X,z.Y,'MS',z.theta);
 assert(max(abs(model.R(:)-z.R(:)))<1e-12);rng(z.expected.rng_before);before=rng;
 x=FDCVV_function(model,[zeros(1,z.b);ones(1,z.b)],z.step,30);after=rng;
 t=load(sprintf('trace_%02d.mat',z.step));tr=t.trace;
 assert(isequal(before,after));assert(isequal(x,z.expected.x));
 assert(isequal(tr.regions,z.expected.regions)&&isequal(tr.K,z.expected.K));
 assert(abs(tr.score-z.expected.score)<1e-12&&max(abs(tr.fast_scores-z.expected.fast_scores))<1e-12);
 assert(isequal(tr.region_winners,z.expected.region_winners));
 expectedK=4-double((z.step-1)/30>=1/3)-double((z.step-1)/30>=1/2);
 assert(tr.K==expectedK&&abs(tr.lambda-(.35+.65*(z.step/30)^(log(.15/.65)/log(.1))))<1e-12);
 rows(end+1,:)=[z.b z.rep z.step tr.K tr.lambda max(abs(x-z.expected.x))];
 fprintf('PASS benchmark=%d rep=%d step=%d K=%d\n',z.b,z.rep,z.step,tr.K);
end
writetable(array2table(rows,'VariableNames',{'benchmark_index','source_repeat','step','K','lambda','max_point_difference'}),fullfile(root,'validation','replay_checks.csv'));
report=struct('status','PASS','cases',size(rows,1),'matlab',version,'toolboxes',{ver},'new_optimized_fits',0,'new_truth_evaluations',0);
fid=fopen(fullfile(root,'validation','FDCVV_replay.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
end
