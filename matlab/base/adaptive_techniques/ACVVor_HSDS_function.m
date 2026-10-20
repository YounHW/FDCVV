function x_new = ACVVor_HSDS_function(obj, A)
% HSDS-ACVVor
% Hierarchical Screening + Directional Search ACVVor
%
% 目标：
% 1) 保留原 ACVVor 的基本逻辑：
%    - 用 LOO + ACE 评价采样价值
%    - 用 Voronoi 单元进行局部搜索
%    - 保留最小距离约束和逐步放宽机制
%
% 2) 减少计算量：
%    - 先对所有 Voronoi 单元做粗筛
%    - 仅对 top-K 单元做 ACE 精筛
%    - 在最佳单元中用定向搜索代替大量随机候选点
%
% 输入：
% obj - Ordinary Kriging 模型对象
% A   - 设计空间范围，A(1,:) 为下界，A(2,:) 为上界
%
% 输出：
% x_new - 新的采样点
%
% 依赖：
% - OK_model 类
% - OK_model/predict 方法
%
% 说明：
% - 这里为保证"可直接运行”，内部自带了 Voronoi 粗支持点生成函数 
% - 若你已有 randomVoronoi，也可替换 buildCoarseVoronoi()

%% ---------------------------
% 设计空间上下界
% ----------------------------
lb = A(1,:);
ub = A(2,:);

m = obj.m;
dim = size(obj.X,2);

%% ---------------------------
% 参数设置
% ----------------------------
lambda_c = 0.7;              % 粗筛：LOO vs 单元稀疏度
lambda   = 0.8;              % 精筛：LOO vs ACE代表点均值
mu       = 0.2;              % 点级评分：ACE vs 距中心距离

K = max(3, ceil(0.10 * m));  % 进入精筛的单元数
q = 5;                       % 每个热点单元代表点数
s_max = 6;                   % 最多搜索方向数
eta_list = [0.25, 0.50, 0.75, 1.00];  % 各方向试探步长比例

% 粗粒度 Voronoi 支持点总数
nSupport = max(200, 30 * m);

%% ---------------------------
% 构造粗粒度 Voronoi 单元
% C{i,1} = 第 i 个单元中心样本点
% C{i,2} = 分配到第 i 个单元的粗支持点
% ----------------------------
C = buildCoarseVoronoi(obj.X, lb, ub, nSupport);

%% ---------------------------
% 基础量
% ----------------------------
% 保持原算法风格：使用原始 LOO（删一点重训）
% 若后续你想进一步加速，可替换成 eLOO_fast_function(obj)
eLOO = eL00_function(obj);

[alpha, d] = ACE_alpha(obj);

%% ==========================================================
% Step 1: 全单元粗筛
% 粗指标：
% 1) eLOO(i)
% 2) Voronoi 单元外扩半径（基于粗支持点）
% ==========================================================
cellRadius = zeros(m,1);
cellNonEmpty = false(m,1);

for i = 1:m
    Xi = C{i,2};
    Pc = C{i,1};
    
    if ~isempty(Xi)
        cellNonEmpty(i) = true;
        dist_center = rowDistances(Xi, Pc);
        cellRadius(i) = max(dist_center);
    else
        cellRadius(i) = 0.0;
    end
end

eLOO_n = normalizeVector(eLOO);
cellRadius_n = normalizeVector(cellRadius);

coarseScore = lambda_c * eLOO_n + (1-lambda_c) * cellRadius_n;
coarseScore(~cellNonEmpty) = -inf;

[~, p_coarse] = sort(coarseScore, 'descend');

nNonEmpty = sum(cellNonEmpty);
if nNonEmpty == 0
    % 极端兜底：如果没有任何非空单元，则直接全局随机选一点
    x_new = lb + rand(1,dim) .* (ub - lb);
    return
end

K = min(K, nNonEmpty);
topCells = p_coarse(1:K);

%% ==========================================================
% Step 2: 仅对 top-K 单元做 ACE 精筛
% 每个热点单元只选 q 个代表点做 ACE 均值近似
% ==========================================================
cellACE_rep = -inf(m,1);
repPoints = cell(m,1);
repACE = cell(m,1);

for ii = 1:numel(topCells)
    i = topCells(ii);
    
    Xi = C{i,2};
    Pc = C{i,1};
    
    Xrep = chooseRepresentativePoints(Xi, Pc, q);
    repPoints{i} = Xrep;
    
    if isempty(Xrep)
        cellACE_rep(i) = -inf;
        repACE{i} = [];
        continue
    end
    
    Erep = zeros(size(Xrep,1),1);
    for j = 1:size(Xrep,1)
        Erep(j) = ACVVor_errorField(obj, eLOO, alpha, Xrep(j,:));
    end
    
    repACE{i} = Erep;
    cellACE_rep(i) = mean(Erep);
end

% 只对 top-K 精筛单元做归一化
e_top = eLOO(topCells);
a_top = cellACE_rep(topCells);

e_top_n = normalizeVector(e_top);
a_top_n = normalizeVector(a_top);

fineScore_top = lambda * e_top_n + (1-lambda) * a_top_n;

[~, idx_sort_top] = sort(fineScore_top, 'descend');
refinedCells = topCells(idx_sort_top);

%% ==========================================================
% Step 3: 按精筛结果，依次在单元内做定向搜索
% 若无满足距离约束的点，则逐步放宽阈值
% ==========================================================
d_current = d;
found = false;

while ~found
    
    for iter = 1:numel(refinedCells)
        k = refinedCells(iter);
        
        Pc = C{k,1};
        Xrep = repPoints{k};
        Erep = repACE{k};
        
        if isempty(Xrep)
            continue
        end
        
        % 构造搜索方向
        U = buildDirectionalSet(Pc, Xrep, Erep, obj.X, s_max);
        
        % 沿方向生成少量定向候选点
        Xcand = directionalCandidates(Pc, U, eta_list, lb, ub, obj.X, k);
        
        % 把代表点也加入候选集，提高稳健性
        Xcand = [Xcand; Xrep];
        Xcand = unique(Xcand, 'rows', 'stable');
        
        if isempty(Xcand)
            continue
        end
        
        % 点级评分：保持原 ACVVor 风格
        errCand = zeros(size(Xcand,1),1);
        distCand = zeros(size(Xcand,1),1);
        minDistCand = zeros(size(Xcand,1),1);
        
        for i = 1:size(Xcand,1)
            errCand(i) = ACVVor_errorField(obj, eLOO, alpha, Xcand(i,:));
            distCand(i) = norm(Xcand(i,:) - Pc);
            
            dist_all = rowDistances(obj.X, Xcand(i,:));
            minDistCand(i) = min(dist_all);
        end
        
        errCand_n = normalizeVector(errCand);
        distCand_n = normalizeVector(distCand);
        
        pointScore = mu * errCand_n + (1-mu) * distCand_n;
        
        feasible = minDistCand >= d_current;
        
        if any(feasible)
            feasible_idx = find(feasible);
            [~, idx_local] = max(pointScore(feasible_idx));
            x_new = Xcand(feasible_idx(idx_local), :);
            found = true;
            break
        end
    end
    
    % 若当前阈值下无解，则放宽距离约束
    if ~found
        d_current = d_current / 2;
        
        % 兜底策略
        if d_current < 1e-6
            x_new = fallbackBestCandidate(refinedCells, C, repPoints, repACE, ...
                                          obj, eLOO, alpha, mu, ...
                                          eta_list, lb, ub, s_max);
            found = true;
        end
    end
end

end


%% ========================================================================
%  辅助函数：构造粗粒度 Voronoi 单元
% ========================================================================
function C = buildCoarseVoronoi(X, lb, ub, nSupport)
% 用一批全局随机支持点近似构造 Voronoi 单元
%
% 输出：
% C{i,1} - 第 i 个单元中心样本点 X(i,:)
% C{i,2} - 属于该单元的粗支持点

m = size(X,1);
dim = size(X,2);

% 全局均匀随机支持点
Z = repmat(lb, nSupport, 1) + rand(nSupport, dim) .* repmat((ub-lb), nSupport, 1);

C = cell(m,2);
for i = 1:m
    C{i,1} = X(i,:);
    C{i,2} = zeros(0, dim);
end

% 按最近样本点分配到对应 Voronoi 单元
for i = 1:nSupport
    z = Z(i,:);
    dist_all = rowDistances(X, z);
    [~, idx_min] = min(dist_all);
    C{idx_min,2} = [C{idx_min,2}; z];
end
end


%% ========================================================================
%  辅助函数：代表点选择
% ========================================================================
function Xrep = chooseRepresentativePoints(Xi, Pc, q)
% 从单元粗支持点 Xi 中选 q 个“远且分散”的代表点

if isempty(Xi)
    Xrep = [];
    return
end

n = size(Xi,1);

if n <= q
    Xrep = Xi;
    return
end

dist_center = rowDistances(Xi, Pc);
[~, idx1] = max(dist_center);

selected = idx1;
remaining = setdiff(1:n, selected);

while numel(selected) < q && ~isempty(remaining)
    best_idx = remaining(1);
    best_score = -inf;
    
    for rr = 1:numel(remaining)
        r = remaining(rr);
        
        d_to_sel = zeros(numel(selected),1);
        for j = 1:numel(selected)
            d_to_sel(j) = norm(Xi(r,:) - Xi(selected(j),:));
        end
        
        score = 0.5 * norm(Xi(r,:) - Pc) + 0.5 * min(d_to_sel);
        
        if score > best_score
            best_score = score;
            best_idx = r;
        end
    end
    
    selected = [selected, best_idx];
    remaining = setdiff(remaining, best_idx);
end

Xrep = Xi(selected,:);
end


%% ========================================================================
%  辅助函数：构造定向搜索方向
% ========================================================================
function U = buildDirectionalSet(Pc, Xrep, Erep, Xall, s_max)
% 方向来源：
% 1) 最远代表点方向
% 2) ACE 最大代表点方向
% 3) 其余代表点方向
% 4) 邻近样本排斥方向

dim = size(Pc,2);
U = zeros(0, dim);

if isempty(Xrep)
    return
end

% 1) 最远代表点方向
dist_rep = rowDistances(Xrep, Pc);
[~, idx_far] = max(dist_rep);
u_far = normalizeDirection(Xrep(idx_far,:) - Pc);
if norm(u_far) > 0
    U = [U; u_far];
end

% 2) ACE 最大代表点方向
if ~isempty(Erep)
    [~, idx_ace] = max(Erep);
    u_ace = normalizeDirection(Xrep(idx_ace,:) - Pc);
    if norm(u_ace) > 0
        U = [U; u_ace];
    end
end

% 3) 其余代表点方向
for i = 1:size(Xrep,1)
    u = normalizeDirection(Xrep(i,:) - Pc);
    if norm(u) > 0
        U = [U; u];
    end
end

% 4) 邻近样本排斥方向
d_all = rowDistances(Xall, Pc);
[~, idx_sort] = sort(d_all, 'ascend');

% 去掉自身
idx_sort(d_all(idx_sort) < 1e-12) = [];

k_near = min(3, numel(idx_sort));
if k_near >= 1
    u_repulse = zeros(1, dim);
    for j = 1:k_near
        v = Pc - Xall(idx_sort(j),:);
        nv = norm(v);
        if nv > 0
            u_repulse = u_repulse + v / nv;
        end
    end
    u_repulse = normalizeDirection(u_repulse);
    if norm(u_repulse) > 0
        U = [U; u_repulse];
    end
end

% 去重（近似）
U = unique(round(U, 10), 'rows', 'stable');

% 截断方向数
if size(U,1) > s_max
    U = U(1:s_max, :);
end
end


%% ========================================================================
%  辅助函数：沿方向生成候选点
% ========================================================================
function Xcand = directionalCandidates(Pc, U, eta_list, lb, ub, Xall, cell_id)
% 从中心点 Pc 沿方向 U 生成少量候选点
% 要求：
% 1) 在设计空间内
% 2) 属于指定 Voronoi 单元

dim = size(Pc,2);
Xcand = zeros(0, dim);

if isempty(U)
    return
end

for i = 1:size(U,1)
    u = U(i,:);
    
    tmax = maxStepInBox(Pc, u, lb, ub);
    if tmax <= 0
        continue
    end
    
    for j = 1:numel(eta_list)
        t = eta_list(j) * tmax;
        
        % 若当前步长点不在目标 Voronoi 单元内，则逐步缩小
        accepted = false;
        for shrinkIter = 1:6
            x = Pc + t * u;
            
            if isInsideBox(x, lb, ub) && isInVoronoiCell(x, Xall, cell_id)
                Xcand = [Xcand; x];
                accepted = true;
                break
            end
            
            t = 0.5 * t;
            if t < 1e-12
                break
            end
        end
        
        if ~accepted
            % do nothing
        end
    end
end

if ~isempty(Xcand)
    Xcand = unique(Xcand, 'rows', 'stable');
end
end


%% ========================================================================
%  辅助函数：兜底策略
% ========================================================================
function x_best = fallbackBestCandidate(refinedCells, C, repPoints, repACE, ...
                                        obj, eLOO, alpha, mu, ...
                                        eta_list, lb, ub, s_max)

x_best = [];
best_score = -inf;

for iter = 1:numel(refinedCells)
    k = refinedCells(iter);
    
    Pc = C{k,1};
    Xrep = repPoints{k};
    Erep = repACE{k};
    
    if isempty(Xrep)
        continue
    end
    
    U = buildDirectionalSet(Pc, Xrep, Erep, obj.X, s_max);
    Xcand = directionalCandidates(Pc, U, eta_list, lb, ub, obj.X, k);
    Xcand = [Xcand; Xrep];
    Xcand = unique(Xcand, 'rows', 'stable');
    
    if isempty(Xcand)
        continue
    end
    
    errCand = zeros(size(Xcand,1),1);
    distCand = zeros(size(Xcand,1),1);
    
    for i = 1:size(Xcand,1)
        errCand(i) = ACVVor_errorField(obj, eLOO, alpha, Xcand(i,:));
        distCand(i) = norm(Xcand(i,:) - Pc);
    end
    
    errCand_n = normalizeVector(errCand);
    distCand_n = normalizeVector(distCand);
    pointScore = mu * errCand_n + (1-mu) * distCand_n;
    
    [score_local, idx_local] = max(pointScore);
    if score_local > best_score
        best_score = score_local;
        x_best = Xcand(idx_local,:);
    end
end

if isempty(x_best)
    % 最终兜底：在 refinedCells 中找一个代表点
    for iter = 1:numel(refinedCells)
        k = refinedCells(iter);
        Xrep = repPoints{k};
        if ~isempty(Xrep)
            x_best = Xrep(1,:);
            return
        end
    end
    
    % 再极端兜底：随机点
    dim = size(obj.X,2);
    x_best = lb + rand(1,dim) .* (ub - lb);
end
end


%% ========================================================================
%  辅助函数：判断点是否属于某 Voronoi 单元
% ========================================================================
function flag = isInVoronoiCell(x, Xall, cell_id)
dist_all = rowDistances(Xall, x);
[~, idx_min] = min(dist_all);
flag = (idx_min == cell_id);
end


%% ========================================================================
%  辅助函数：求方向在盒约束中的最大步长
% ========================================================================
function tmax = maxStepInBox(x0, u, lb, ub)
t_list = inf(1, numel(u));

for k = 1:numel(u)
    if abs(u(k)) < 1e-14
        t_list(k) = inf;
    elseif u(k) > 0
        t_list(k) = (ub(k) - x0(k)) / u(k);
    else
        t_list(k) = (lb(k) - x0(k)) / u(k);
    end
end

tmax = min(t_list);

if isempty(tmax) || ~isfinite(tmax) || tmax < 0
    tmax = 0;
end
end


%% ========================================================================
%  工具函数：盒约束检查
% ========================================================================
function flag = isInsideBox(x, lb, ub)
flag = all(x >= lb) && all(x <= ub);
end


%% ========================================================================
%  工具函数：方向归一化
% ========================================================================
function u = normalizeDirection(v)
nv = norm(v);
if nv < 1e-14
    u = zeros(size(v));
else
    u = v / nv;
end
end


%% ========================================================================
%  工具函数：矩阵每行到某点的欧氏距离
% ========================================================================
function d = rowDistances(X, x)
diffMat = X - repmat(x, size(X,1), 1);
d = sqrt(sum(diffMat.^2, 2));
end


%% ========================================================================
%  原 ACVVor 所需辅助函数：LOO 误差
% ========================================================================
function eLOO_perPoint = eL00_function(obj)
% 与你原代码保持一致：删除一点后重新训练 OK 模型

eLOO_perPoint = zeros(obj.m,1);

for i = 1:obj.m
    
    Xp = obj.X;
    Xp(i,:) = [];
    
    Yp = obj.Y;
    Yp(i) = [];
    
    M_hatp = OK_model(obj.auto_correlation_function, Xp, Yp, obj.theta_opti_technique);
    
    y = obj.Y(i);
    [mu_hatp, ~] = predict(M_hatp, obj.X(i,:));
    
    eLOO_perPoint(i) = norm(y - mu_hatp);
end
end


%% ========================================================================
%  原 ACVVor 所需辅助函数：ACE 型累计误差场
% ========================================================================
function E = ACVVor_errorField(obj, eLOO, alpha, x)
E = 0.0;

for i = 1:obj.m
    E = E + eLOO(i) * exp(-alpha * norm(obj.X(i,:) - x));
end
end


%% ========================================================================
%  原 ACVVor 所需辅助函数：ACE 参数
% ========================================================================
function [alpha, d] = ACE_alpha(obj)
min_distances = zeros(obj.m,1);

for i = 1:obj.m
    
    iter = 1;
    distance_min = zeros(obj.m-1,1);
    
    for j = 1:obj.m
        if ~(i == j)
            distance_min(iter) = norm(obj.X(i,:) - obj.X(j,:));
            iter = iter + 1;
        end
    end
    
    min_distances(i) = min(distance_min);
end

max_distance = max(min_distances);
DOI_distance = 0.5 * max_distance;
worst_case_DOI = 1e-5;

if DOI_distance < eps
    alpha = 1.0;
else
    alpha = -log(worst_case_DOI) / DOI_distance;
end

d = 0.5 * mean(min_distances);

if d < eps
    d = 1e-6;
end
end


%% ========================================================================
%  原 ACVVor 所需辅助函数：归一化
% ========================================================================
function v_n = normalizeVector(v)
v = v(:);

finiteMask = isfinite(v);

if ~any(finiteMask)
    v_n = zeros(size(v));
    return
end

vf = v(finiteMask);
vmin = min(vf);
vmax = max(vf);

v_n = zeros(size(v));

if abs(vmax - vmin) < eps
    v_n(finiteMask) = 0;
else
    v_n(finiteMask) = (vf - vmin) ./ (vmax - vmin);
end

% 非有限值统一置 0
v_n(~finiteMask) = 0;
end