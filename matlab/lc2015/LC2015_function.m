function xnew = LC2015_function(model,A,iteration)
% LC2015 s=1 adjusted variance, all Voronoi regions compete.
% Use original project's LHS generator and n*d*2000 candidate density.
candidate_rng=rng;
[xnew,info,search]=lc_search(model,A);
trace=struct('iteration',iteration,'parent_n',model.m,'theta',model.theta, ...
    'candidate_rng',candidate_rng,'candidate_count',search.candidate_count, ...
    'selected_point',xnew,'selected_index',info.index, ...
    'selected_owner',info.owner(info.index),'selected_score',info.score(info.index), ...
    'loo',info.loo,'ratio',info.ratio,'search',search);
% Candidate pool and full scores can be exactly regenerated from saved RNG.
save(sprintf('selection_%02d.mat',iteration),'trace');
end

