function x_new = FDCVV_function(obj,A,m,mMax)
% Adopted FDCVV. K=4/3/2 at completed progress 1/3 and 1/2.
% lambda=.35+.65*p^alpha, alpha=log(.15/.65)/log(.10), p=m/mMax.
% Within-cell variance and nearest-sample distance normalization.
% ABL_B is benchmark index 1..3; ABL_REP is source repeat 106..135.
% Candidate stream is restored on return. trace_XX.mat records selection.
global ABL_B ABL_REP; prior_rng=rng; restore_rng=onCleanup(@()rng(prior_rng)); candidate_seed=271010000+ABL_B*100000+ABL_REP*100+m; rng(candidate_seed,'twister'); 
lb = A(1,:);
ub = A(2,:);

% ===================== ablation selector =====================
% 可选：'Full','NoLOOCVSeeds','NoVoronoi','FixedLambda','NoDynamicTopK','OnlySigmaInScore','OnlyDistInScore'
ablation_mode = 'Full';
% ============================================================

% -------------------- fixed defaults (Full) --------------------
useLOOCVSeeds   = true;   % M1 
useVoronoi      = true;   % M2

lambdaMode     = 'dynamic'; % M3: 'dynamic'|'fixed'
lambdaFixed    = 0.5;

topKMode        = 'dynamic'; % M-topK: 'dynamic'|'fixed'
kMax           = 3;
p_target       = 0.5;
fixedTopK      = 2;

useSigmaInScore = true;  % M4
useDistInScore  = true; % M4

doNormalization  = true; % M5
epsVal = 1e-12;
% nearestX distance term (fixed): dist(x)=min distance to existing samples

% -------------------- apply ablation --------------------
switch ablation_mode
    case 'Full'
        % keep defaults
    case 'NoLOOCVSeeds'
        useLOOCVSeeds = false; % seed selection random
    case 'NoVoronoi'
        useVoronoi = false;    % strong: no Voronoi cells & no Voronoi-based candidate partition
    case 'FixedLambda'
        lambdaMode = 'fixed';
        lambdaFixed = 0.5;
    case 'NoDynamicTopK'
        topKMode = 'fixed';
        fixedTopK = 2;
    case 'OnlySigmaInScore'
        useSigmaInScore = true;
        useDistInScore  = false;
    case 'OnlyDistInScore'
        useSigmaInScore = false;
        useDistInScore  = true;
    otherwise
        error('Unknown ablation_mode: %s', ablation_mode);
end

% -------------------- dynamic lambda (M3) --------------------
lambda0 = 0.35;
p = m / max(1, mMax);
p = min(max(p,0),1);  % clamp to [0,1]
p_cross = 0.10;

alpha = log((0.5 - lambda0) / (1 - lambda0)) / log(p_cross);

switch lower(lambdaMode)
    case 'dynamic'
        lambda = lambda0 + (1 - lambda0) * (p^alpha);
    case 'fixed'
        lambda = lambdaFixed;
    otherwise
        error('Unknown lambdaMode');
end

% =================== strong NoVoronoi ===================
if ~useVoronoi
    % --- candidate generation: 全域 LHS 候选，不构造 Voronoi cells ---
    n_of_existing_P = size(obj.X,1);
    dimension = size(obj.X,2);
    w = 1000;

    nCand = n_of_existing_P * dimension * w;
    if nCand <= 0
        x_new = lb;
        return;
    end

    xn = lhs_scaled(nCand, lb, ub);

    % 原 randomVoronoi 的全局候选等价于：中心点(obj.X) + 全部xn
    % 这里按你的要求做“去重（unique rows）”
    globalCand = [obj.X; xn];
    globalCand = unique(globalCand, 'rows');

    if isempty(globalCand)
        x_new = lb;
        return;
    end

    % --- scoring on all candidates (global normalization) ---
    [~, sigma2] = obj.predict(globalCand);
    sigma2 = sigma2(:);

    % nearestX distance term (fixed)
    distToCenter = min(pdist2(globalCand, obj.X), [], 2);

    % normalization (per global batch)
    if doNormalization
        sigma2_n = (sigma2 - min(sigma2)) / (max(sigma2) - min(sigma2) + epsVal);
        dist_n   = (distToCenter - min(distToCenter)) / (max(distToCenter) - min(distToCenter) + epsVal);
    else
        sigma2_n = sigma2;
        dist_n   = distToCenter;
    end

    % score composition
    if useSigmaInScore && useDistInScore
        score = lambda * sigma2_n + (1 - lambda) * dist_n;
    elseif useSigmaInScore
        score = sigma2_n;
    elseif useDistInScore
        score = dist_n;
    else
        error('Invalid score flags: both are false.');
    end

    [~, idx] = max(score);
    x_new = globalCand(idx,:);
    trace=struct('completed_progress',(m-1)/mMax,'step',m,'budget',mMax,'p',p,'lambda',lambda,'regions',[],'K',NaN,'x',x_new,'score',score(idx),'mode',ablation_mode);
    if exist('region_winners','var'),trace.region_winners=region_winners;end;trace.candidate_seed=candidate_seed;trace.rng_before=prior_rng;save(sprintf('trace_%02d.mat',m),'trace');
    return;
end

% =================== Voronoi mode (original logic) ===================
% -------------------- Voronoi candidates --------------------
C  = randomVoronoi(obj.X,lb,ub);
numCells = size(C,1);

% -------------------- LOOCV sensitivity per seed point (optional) --------------------
eLOO_fast_perPoint = [];
if useLOOCVSeeds
    eLOO_fast_perPoint = sqrt(eLOO_fast_function(obj)); % squared -> sqrt
    numSeeds = min(numCells, numel(eLOO_fast_perPoint));
else
    numSeeds = numCells;
end

if numSeeds <= 0
    x_new = lb;
    return;
end

% -------------------- choose topK sensitive/random seeds --------------------
% K uses completed additions; lambda retains its original progress.
pK=(m-1)/mMax;
topK=min(numSeeds,4 - double((m-1)/mMax >= 1/3) - double((m-1)/mMax >= 1/2));
[~,idxSorted]=sort(eLOO_fast_perPoint(1:numSeeds),'descend');
topSeedIdx=idxSorted(1:topK);

% -------------------- candidate selection & scoring --------------------
region_winners=[]; bestScore = -inf;
x_new = [];

for t = 1:numel(topSeedIdx)
    k = topSeedIdx(t); % cell index
    if k < 1 || k > numCells, continue; end

    cellCand = C{k,2};
    if isempty(cellCand), continue; end

    [~, sigma2] = obj.predict(cellCand);
    sigma2 = sigma2(:);

    distToCenter = min(pdist2(cellCand, obj.X), [], 2);

    if doNormalization
        sigma2_n = (sigma2 - min(sigma2)) / (max(sigma2) - min(sigma2) + epsVal);
        dist_n   = (distToCenter - min(distToCenter)) / (max(distToCenter) - min(distToCenter) + epsVal);
    else
        sigma2_n = sigma2;
        dist_n   = distToCenter;
    end

    if useSigmaInScore && useDistInScore
        score = lambda * sigma2_n + (1 - lambda) * dist_n;
    elseif useSigmaInScore
        score = sigma2_n;
    elseif useDistInScore
        score = dist_n;
    else
        error('Invalid score flags: both are false.');
    end

    [cellBestScore, idx] = max(score); region_winners=[region_winners; k,cellCand(idx,:),sigma2(idx),distToCenter(idx),sigma2_n(idx),dist_n(idx),cellBestScore];
    if cellBestScore > bestScore
        bestScore = cellBestScore;
        x_new = cellCand(idx,:);
    end
end

if isempty(x_new)
    x_new = lb;
end
trace=struct('completed_progress',(m-1)/mMax,'step',m,'budget',mMax,'p',p,'lambda',lambda,'regions',topSeedIdx,'K',numel(topSeedIdx),'x',x_new,'score',bestScore,'mode',ablation_mode,'fast_scores',eLOO_fast_perPoint);
if exist('region_winners','var'),trace.region_winners=region_winners;end;trace.candidate_seed=candidate_seed;trace.rng_before=prior_rng;save(sprintf('trace_%02d.mat',m),'trace');
end


% -------------------- LOOCV fast per point --------------------
function eLOO_perPoint = eLOO_fast_function(obj)
% Fast LOOCV per point (squared LOOCV error) for Ordinary Kriging
    beta_hat = compute_beta_hat(obj, obj.R);
    H = obj.F * ((obj.F' * obj.F) \ obj.F');
    d = obj.Y - obj.F * beta_hat;

    inv_R = inv(obj.R);

    eLOO_perPoint = zeros(obj.m,1);
    for i = 1:obj.m
        eLOO_perPoint(i) = ( (inv_R(i,:) * ( d + H(:,i) * (d(i)/(1-H(i,i))) )) ...
                             / (inv_R(i,i)) )^2;
    end
end



