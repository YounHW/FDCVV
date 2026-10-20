function [xnew,info,search] = lc_search(model,A)
% Nested independent LHS batches; retain best candidate across all batches.
% Search convergence concerns the acquisition maximum, not test error.
batch=model.m*size(model.X,2)*2000;best=-inf;stable=0;history=[];
for k=1:16
 C=lhs_scaled(batch,A(1,:),A(2,:));[x,trial]=lc2015_select(model,C);
 value=trial.score(trial.index);old=best;
 if value>best,best=value;xnew=x;info=trial;winningBatch=k;end
 gain=(best-old)/max(abs(old),realmin);
 if k>1 && gain<=0.001,stable=stable+1;else,stable=0;end
 history(k,:)=[k value best gain]; %#ok<AGROW>
 if k>=4 && stable>=3,break;end
end
search=struct('batches',k,'batch_size',batch,'candidate_count',batch*k, ...
 'winning_batch',winningBatch,'history',history,'converged',stable>=3, ...
 'rule','min4 max16 batches; 3 successive relative improvements <=0.1%');
assert(search.converged,'LC search reached maximum batches without convergence');
end
