function x_new = FDCVV_function(obj,A,m,mMax)
% Fast Dynamic Cross-Validation Voronoi (FDCVV)
%
% 功能概述
% --------
% 该函数用在给定决策于变量边界 [lb, ub] 内，从候选集合中选择一个下一步采样点 x_new。
% 选择准则由以下几部分组成（可通过 ablation_mode 控制消融）：
%   1) 预测不确定性（sigma2，来自 obj.predict）
%   2) 候选点与已有样本的“空间差异/距离”（distToCenter 或到最近已有点距离）
%   3) dynamic λ：权重随当前迭代阶段（由 m, mMax 推断）自适应调整
%   4) topK：在 Voronoi cell 中选择 topK 个“最敏感/最可能更优”的 cell 进行精评分支
%   5) （可选）LOOCV seeds：用快速 LOOCV 指标对每个 Voronoi seed 点进行敏感性排序
%
% 在 **Voronoi mode** 下，会先将空间划分为 Voronoi cells（由 randomVoronoi(obj.X,lb,ub) 生成）。
% 然后在选中的 topK cells 内，对该 cell 的候选点集 cellCand 逐点预测并打分，取分数最大者作为 x_new。
%
% 在 **NoVoronoi 模式**（useVoronoi=false）下，不构造 Voronoi cells，而是对“全域候选”
% （全域 LHS 采样生成的候选，与 obj.X 合并并 unique 去重后）直接做全局打分取最大者。
%
% 用途
% ----
% 适用于（通常）基于高斯过程/普通克里金（Ordinary Kriging）框架的主动学习/自适应采样：
% 以预测方差 sigma2 与空间覆盖性距离为目标，并辅以 LOOCV 指导的候选筛选来提高效率并降低误差风险。
%
% 输入参数
% ----------
% obj : 结构体/对象（依赖字段与方法）
%   必须包含（至少）以下成员/方法：
%     - obj.X : 已有样本点，大小为 [m0, d]
%               其中 d 为输入维度，m0 为已有样本数
%     - obj.predict(Z) : 对输入候选点 Z（大小 [n, d]）返回预测结果
%               其中 sigma2 是预测方差（代码中使用第二个输出）
%               要求： [mu, sigma2] = obj.predict(Z)
%               且 sigma2 可 reshape 成列向量 [n,1]
%     - obj.m : 代码内 eLOO_fast_function 使用的训练样本数（用于循环次数）
%               理想情况下 obj.m == size(obj.X,1)
%     - obj.R : 克里金相关矩阵（通常为 [obj.m, obj.m]）
%               要求：可逆（代码中使用 inv(obj.R)）
%     - obj.F : 设计矩阵（通常为常基/回归基函数），用于计算 beta_hat
%               维度需满足代码计算：
%                 beta_hat = compute_beta_hat(obj, obj.R);
%                 H = obj.F * ((obj.F' * obj.F) \ obj.F');
%                 以及 d = obj.Y - obj.F * beta_hat;
%     - obj.Y : 观测响应向量，维度为 [obj.m, 1]
%
% A : 2×d 的边界矩阵
%   A(1,:) 为下界 lb
%   A(2,:) 为上界 ub
%
% m : 当前阶段的样本数/迭代计数（用于 dynamic λ）
%   代码中使用 p = m / max(1,mMax)，并将其 clamp 到 [0,1]
%
% mMax : 允许的最大样本数/最大迭代数
%
% 输出参数
% ----------
% x_new : 1×d 的下一步建议采样点（取分数 score 最大对应的候选点）
%         若出现候选为空或异常路径，则退化返回 lb
%
% 关键内部设置（可通过 ablation_mode 影响）
% --------------------------------------
% ablation_mode = 'Full' (默认)
% 允许的消融选项：
%   - 'Full'
%   - 'NoLOOCVSeeds'     : 关闭 LOOCV seeds，用随机 seeds 或直接 topK
%   - 'NoVoronoi'         : 完全禁用 Voronoi cells，改用全域候选直接全局评分
%   - 'FixedLambda'      : λ 改为固定 lambdaFixed
%   - 'NoDynamicTopK'    : topK 改为固定 fixedTopK
%   - 'OnlySigmaInScore' : score 只使用 sigma2 项（不使用距离）
%   - 'OnlyDistInScore'  : score 只使用距离项（不使用 sigma2）
%
% dynamic lambda (M3)
% -----------------
% 使用参数：
%   lambda0 = 0.25
%   p = m / max(1,mMax)  -> clamp 到 [0,1]
%   p_cross = 0.20
% alpha 通过满足一定交叉点条件推导：
%   alpha = log((0.5 - lambda0) / (1 - lambda0)) / log(p_cross)
% 最终：
%   lambda = lambda0 + (1 - lambda0) * (p^alpha)     (dynamic 模式)
% 或 lambda = lambdaFixed (fixed 模式)
%
% topK 动态选择 (M-topK)
% ---------------------
% kMax = 3
% p_target = 0.5
% 使用 kMaxEff = min(kMax, numSeeds)
% 当 kMaxEff > 1：
%   gamma = ceil( log(1/(kMaxEff-1)) / log(1 - p_target) )
%   term  = (kMaxEff - 1) * ((1 - p)^gamma)
%   topK  = 1 + floor(term)，并夹到 [1, kMaxEff]
% 然后：
%   若使用 LOOCV seeds：按 eLOO_fast_perPoint 降序选择 topK
%   否则随机选择 topK
%
% 评分 score 的定义
% ----------------
% 对每个候选点，计算：
%   - sigma2：预测方差（obj.predict 输出）
%   - distToCenter：到已有样本的最近距离（用 min(pdist2(候选, obj.X))）
%
% 归一化（doNormalization=true）
%   sigma2_n = (sigma2 - min(sigma2)) / (max(sigma2) - min(sigma2) + epsVal)
%   dist_n   = (distToCenter - min(distToCenter)) / (max(distToCenter) - min(distToCenter) + epsVal)
%
% 组合：
%   若 useSigmaInScore && useDistInScore：
%       score = lambda * sigma2_n + (1 - lambda) * dist_n
%   仅 sigma 或仅 dist 则相应简化
%
% Voronoi 分支（useVoronoi=true）
% ----------------------------
% 1) C = randomVoronoi(obj.X,lb,ub)
%    期望：C 为 cell 数组，且 C{k,2} 为第 k 个 cell 的候选集合 cellCand
% 2) 可选 LOOCV：
%    eLOO_fast_perPoint = sqrt(eLOO_fast_function(obj))
%    numSeeds = min(numCells, numel(eLOO_fast_perPoint))
% 3) 选 topK 个最敏感 seeds 对应的 Voronoi cells
% 4) 对每个 selected cellCand：
%    - 预测 sigma2
%    - 计算 distToCenter
%    - 计算 score
%    - 取 cell 内最大分数点作为该 cell 的最佳点
% 5) 全局比较各 cell 最佳点，选最终分数最大者 x_new
%
% NoVoronoi 分支（useVoronoi=false）
% --------------------------------
% 1) 用全域 LHS 生成候选：
%    nCand = n_of_existing_P * dimension * w, w=1000
%    xn = lhs_scaled(nCand, lb, ub)
% 2) 与已有样本合并并去重：
%    globalCand = unique([obj.X; xn],'rows')
% 3) 对 globalCand 全量预测 sigma2 并计算最近距离 distToCenter
% 4) 全量归一化并按 score 加权取最大者
%
% LOOCV 快速计算 eLOO_fast_function
% -------------------------------
% eLOO_fast_function(obj) 计算普通克里金（Ordinary Kriging）：
%   每个样本点 i 的 squared LOOCV error（平方的 LOOCV 误差）
% 其依赖：
%   - beta_hat = compute_beta_hat(obj, obj.R)
%   - H = obj.F * ((obj.F' * obj.F) \ obj.F')
%   - d = obj.Y - obj.F * beta_hat
%   - inv_R = inv(obj.R)
% 然后使用类似于 Kriging 的 hat-matrix/解析式更新，避免逐点重拟合。
% 最外层代码再对其开方 sqrt(eLOO_fast_function(obj)) 以得到“幅值”形式 eLOO_fast_perPoint。
%
% 注意事项
% --------
% 1) 代码中使用 inv(obj.R)，数值稳定性与性能可能受 obj.R 条件数影响。
%    若 obj.R 近似奇异，可能导致 NaN/Inf 或结果不稳定。
% 2) Voronoi mode 依赖 randomVoronoi 的返回格式：
%    必须保证 C{k,2} 非空且为 [nCellCand, d] 的候选点矩阵。
% 3) NoVoronoi 模式会生成很大量候选（w=1000），predict 的计算成本可能较高。
%
% 示例（伪代码）
% --------------
%   lb = A(1,:);
%   ub = A(2,:);
%   x_next = FDCVV_function(obj,A,m,mMax);
%
% 作者/版本
% --------
%   根据你的原始实现整理注释（请按需要补充作者信息与版本号）
%

% ====== 下面是你原始代码开始 ======

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
lambda0 = 0.25;
p = m / max(1, mMax);
p = min(max(p,0),1);  % clamp to [0,1]
p_cross = 0.20;

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
if strcmpi(topKMode,'dynamic')
    kMaxEff = min(kMax, numSeeds);
    if kMaxEff <= 1
        topSeedIdx = 1;
    else
        gamma = ceil( log(1/(kMaxEff-1)) / log(1 - p_target) );
        term  = (kMaxEff - 1) * ((1 - p) ^ gamma);
        topK   = 1 + floor(term);
        topK   = max(1, min(topK, kMaxEff));

        if useLOOCVSeeds
            [~, idxSorted] = sort(eLOO_fast_perPoint(1:numSeeds), 'descend');
            topSeedIdx = idxSorted(1:topK);
        else
            topSeedIdx = randperm(numSeeds, topK);
        end
    end
else
    topK = min(max(1, fixedTopK), numSeeds);
    if useLOOCVSeeds
        [~, idxSorted] = sort(eLOO_fast_perPoint(1:numSeeds), 'descend');
        topSeedIdx = idxSorted(1:topK);
    else
        topSeedIdx = randperm(numSeeds, topK);
    end
end

% -------------------- candidate selection & scoring --------------------
bestScore = -inf;
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

    [cellBestScore, idx] = max(score);
    if cellBestScore > bestScore
        bestScore = cellBestScore;
        x_new = cellCand(idx,:);
    end
end

if isempty(x_new)
    x_new = lb;
end
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
