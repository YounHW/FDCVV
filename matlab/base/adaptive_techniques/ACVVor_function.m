%% ACE-guided CVVor (ACVVor)
% 融合式自适应采样算法：
% 1) ACE 的“累计误差场”思想
% 2) CVVor 的“Voronoi 单元局部填充”思想
%
% 核心思路：
% - 先利用 LOO 误差 + ACE 误差场，对每个 Voronoi 单元进行评分
% - 再在评分最高的 Voronoi 单元中，选择“误差大且离中心较远”的点
%
% 输入：
% obj - Ordinary Kriging 模型对象
% A   - 设计空间范围，A(1,:) 为下界，A(2,:) 为上界
%
% 输出：
% x_new - 新的采样点

function x_new = ACVVor_function(obj,A)

% 添加辅助函数路径


% 提取设计空间上下界
lb = A(1,:);
ub = A(2,:);

% 基于当前已有样本点 obj.X 构造随机 Voronoi 划分
% C 一般是一个 cell 数组：
% C{i,1} - 第 i 个 Voronoi 单元对应的中心样本点
% C{i,2} - 第 i 个 Voronoi 单元中的候选点集合
C = randomVoronoi(obj.X,lb,ub);

% 计算每个已有样本点的留一交叉验证误差（LOO error）
% eLOO(i) 越大，说明第 i 个样本点附近模型越不可靠
eLOO = eL00_function(obj);

% 计算 ACE 中使用的两个关键参数：
% alpha - 误差在空间中扩散/衰减的速度
% d     - 新点与已有样本点之间的最小距离阈值
[alpha,d] = ACE_alpha(obj);

% ---------------------------
% 参数设置：两个融合权重
% ---------------------------

% lambda 用于“Voronoi 单元评分”
% lambda = 1   -> 完全依赖局部 LOO 误差
% lambda = 0   -> 完全依赖单元内 ACE 平均误差
lambda = 0.8;

% mu 用于“单元内候选点评分”
% mu = 1   -> 完全依赖 ACE 误差场
% mu = 0   -> 完全依赖 CVVor 的最远点思想
mu = 0.2;

% 样本点总数
m = obj.m;

% cellACE(i) 用来保存第 i 个 Voronoi 单元的平均 ACE 误差场值
cellACE = zeros(m,1);

% ==========================================================
% 第一步：计算每个 Voronoi 单元的 ACE 平均误差场值
% ==========================================================
for i=1:m
    
    % 取出第 i 个 Voronoi 单元中的候选点集合
    Xi = C{i,2};
    
    % 用于累计该单元内所有候选点的 ACE 误差场值
    val = 0.0;
    
    % 遍历当前 Voronoi 单元中的每一个候选点
    for j=1:size(Xi,1)
        % 计算候选点 Xi(j,:) 处的 ACE 累计误差场值
        val = val + ACVVor_errorField(obj,eLOO,alpha,Xi(j,:));
    end
    
    % 如果该单元非空，则取平均值
    if size(Xi,1) > 0
        cellACE(i) = val/size(Xi,1);
    else
        % 如果为空，则该单元评分记为 0
        cellACE(i) = 0.0;
    end
end

% ==========================================================
% 第二步：对单元评分所用的两个量进行归一化
% ==========================================================
% eLOO_n    - 每个样本点的局部 LOO 误差归一化
% cellACE_n - 每个 Voronoi 单元平均 ACE 值归一化
eLOO_n = normalizeVector(eLOO);
cellACE_n = normalizeVector(cellACE);

% ==========================================================
% 第三步：计算每个 Voronoi 单元的综合评分
% ==========================================================
% 单元评分 = lambda * 局部 LOO 误差 + (1-lambda) * 单元平均 ACE 误差
% 这样既考虑“局部最差点”，也考虑“全局累计误差场”
cellScore = lambda * eLOO_n + (1-lambda) * cellACE_n;

% 按照单元评分从大到小排序
% p(1) 对应当前最值得采样的 Voronoi 单元编号
[~,p] = sort(cellScore,'descend');

% 当前最小距离阈值，初始为 ACE 计算得到的 d
d_current = d;

% found 用于标记是否已找到满足要求的新点
found = false;

% ==========================================================
% 第四步：按 Voronoi 单元评分高低依次尝试选点
% 如果没有满足距离约束的点，则逐步放宽距离阈值
% ==========================================================
while ~found
    
    % 依次检查排序后的每个 Voronoi 单元
    for iter = 1:m
        
        % 当前被选择的 Voronoi 单元编号
        k = p(iter);
        
        % 该单元的中心样本点（CVVor 中的敏感点）
        P_sensitive = C{k,1};
        
        % 当前 Voronoi 单元中的候选点集合
        Xcand = C{k,2};
        
        % 若该单元没有候选点，跳过
        if isempty(Xcand)
            continue
        end
        
        % errCand(i)     - 候选点的 ACE 误差场值
        % distCand(i)    - 候选点到该单元中心点的距离
        % minDistCand(i) - 候选点到所有已有样本点中的最小距离
        errCand = zeros(size(Xcand,1),1);
        distCand = zeros(size(Xcand,1),1);
        minDistCand = zeros(size(Xcand,1),1);
        
        % 遍历当前单元中的每一个候选点
        for i=1:size(Xcand,1)
            
            % 计算该候选点的 ACE 累计误差场值
            errCand(i) = ACVVor_errorField(obj,eLOO,alpha,Xcand(i,:));
            
            % 计算该候选点到当前单元中心样本点的距离
            % 这是借鉴 CVVor 的“最远点优先”思想
            distCand(i) = norm(Xcand(i,:) - P_sensitive);
            
            % 计算该候选点到所有已有样本点的距离
            dist_all = zeros(obj.m,1);
            for j=1:obj.m
                dist_all(j) = norm(Xcand(i,:) - obj.X(j,:));
            end
            
            % 取最小距离，用于判断是否太靠近已有样本点
            minDistCand(i) = min(dist_all);
        end
        
        % 对候选点的两个评分量进行归一化
        errCand_n = normalizeVector(errCand);
        distCand_n = normalizeVector(distCand);
        
        % ------------------------------------------------------
        % 单元内部候选点综合评分
        % pointScore 越大，说明该点越值得被选中
        %
        % 组成：
        % mu * 误差场值 + (1-mu) * 与中心点距离
        %
        % 含义：
        % - 第一项鼓励去累计误差高的地方
        % - 第二项鼓励去该 Voronoi 单元中更“外侧”的位置
        % ------------------------------------------------------
        pointScore = mu * errCand_n + (1-mu) * distCand_n;
        
        % 满足最小距离约束的候选点
        % 即候选点与任意已有样本点的距离都至少为 d_current
        feasible = minDistCand >= d_current;
        
        % 如果存在可行候选点
        if any(feasible)
            feasible_idx = find(feasible);
            
            % 在所有可行点中，选择综合评分最大的那个
            [~,idx_local] = max(pointScore(feasible_idx));
            x_new = Xcand(feasible_idx(idx_local),:);
            
            % 标记找到解并退出
            found = true;
            break
        end
    end
    
    % ------------------------------------------------------
    % 如果当前所有 Voronoi 单元都找不到满足距离约束的候选点
    % 则逐步放宽距离要求
    % ------------------------------------------------------
    if ~found
        d_current = d_current / 2;
        
        % 如果阈值已经很小了，说明很难再严格满足距离约束
        % 则采用兜底策略：在评分最高的单元中，选综合评分最高的点
        if d_current < 1e-6
            
            % 使用评分最高的 Voronoi 单元
            k = p(1);
            P_sensitive = C{k,1};
            Xcand = C{k,2};
            
            errCand = zeros(size(Xcand,1),1);
            distCand = zeros(size(Xcand,1),1);
            
            for i=1:size(Xcand,1)
                % 候选点 ACE 误差场值
                errCand(i) = ACVVor_errorField(obj,eLOO,alpha,Xcand(i,:));
                
                % 候选点到单元中心距离
                distCand(i) = norm(Xcand(i,:) - P_sensitive);
            end
            
            % 归一化
            errCand_n = normalizeVector(errCand);
            distCand_n = normalizeVector(distCand);
            
            % 计算综合评分
            pointScore = mu * errCand_n + (1-mu) * distCand_n;
            
            % 选评分最高的点
            [~,idx_best] = max(pointScore);
            x_new = Xcand(idx_best,:);
            
            found = true;
        end
    end
end

end


function eLOO_perPoint = eL00_function(obj)
% 计算每个样本点的留一交叉验证误差（LOO error）
%
% 对第 i 个样本点：
% 1. 将该点从训练集中删除
% 2. 用剩余样本重新训练一个 Ordinary Kriging 模型
% 3. 用新模型预测被删掉的那个点
% 4. 真实值与预测值之差，即为该点的 LOO 误差
%
% 输入：
% obj - Ordinary kriging 模型对象
%
% 输出：
% eLOO_perPoint - 每个样本点对应的 LOO 误差

eLOO_perPoint = zeros(obj.m,1);

for i=1:obj.m
    
    % 复制训练输入并删除第 i 个样本点
    Xp = obj.X;
    Xp(i,:) = [];
    
    % 复制训练输出并删除第 i 个样本值
    Yp = obj.Y;
    Yp(i) = [];
    
    % 用删点后的数据重新训练一个 Kriging 模型
    M_hatp = OK_model(obj.auto_correlation_function,Xp,Yp,obj.theta_opti_technique);
    
    % 第 i 个样本点真实输出值
    y = obj.Y(i);
    
    % 用删点后的模型在被删掉的点处做预测
    [mu_hatp,~] = predict(M_hatp,obj.X(i,:));
    
    % 计算该点的 LOO 误差
    % 对标量来说，norm 等价于绝对值
    eLOO_perPoint(i) = norm(y - mu_hatp);
end
end


function E = ACVVor_errorField(obj,eLOO,alpha,x)
% ACE 型累计误差场
%
% 给定任意候选点 x，计算它受到所有已有样本点误差的累计影响：
%
% E(x) = sum_i eLOO(i) * exp(-alpha * ||x_i - x||)
%
% 含义：
% - eLOO(i) 越大，说明第 i 个样本点附近误差越大
% - 距离越近，该样本点对 x 的影响越强
% - 距离越远，影响按指数形式衰减
%
% 输入：
% obj   - Ordinary kriging 模型对象
% eLOO  - 所有样本点的 LOO 误差
% alpha - 衰减参数
% x     - 候选点
%
% 输出：
% E - x 处的累计误差场值

E = 0.0;

for i=1:obj.m
    E = E + eLOO(i) * exp(-alpha * norm(obj.X(i,:) - x));
end
end


function [alpha, d] = ACE_alpha(obj)
% 计算 ACE 中使用的两个参数：
%
% alpha - 误差扩散的衰减系数
% d     - 新点与已有样本点的最小距离阈值
%
% 具体做法：
% 1. 对每个样本点，计算其最近邻距离
% 2. 用最近邻距离估计空间尺度
% 3. 根据经验公式得到 alpha 和 d
%
% 输入：
% obj - Ordinary kriging 模型对象
%
% 输出：
% alpha - ACE 误差衰减参数
% d     - 最小距离阈值

% min_distances(i) 表示第 i 个样本点到其最近邻样本点的距离
min_distances= zeros(obj.m,1);

for i=1:obj.m
    
    clear distance_min
    iter = 1;
    
    % 第 i 个点与其他所有点的距离，共 obj.m-1 个
    distance_min = zeros(obj.m-1,1);
    
    for j=1:obj.m
        if ~(i==j)
            distance_min(iter) = norm(obj.X(i,:) - obj.X(j,:));
            iter = iter +1;
        end
    end
    
    % 保存第 i 个点的最近邻距离
    min_distances(i) = min(distance_min);
end

% 所有样本点最近邻距离中的最大值
max_distance = max(min_distances);

% DOI_distance: influence distance，可理解为一个影响半径尺度
DOI_distance = 0.5 * max_distance;

% worst_case_DOI 表示在 DOI_distance 处，误差影响衰减到的很小值
worst_case_DOI = 10^(-5);

% 根据 exp(-alpha * DOI_distance) = worst_case_DOI 反解 alpha
alpha = -log(worst_case_DOI) / DOI_distance;

% d 取最近邻距离均值的一半，作为新点与旧点的最小安全距离
d = 0.5 * mean(min_distances);

end


function v_n = normalizeVector(v)
% 将向量归一化到 [0,1]
%
% 输入：
% v - 待归一化向量
%
% 输出：
% v_n - 归一化后的向量
%
% 如果向量所有元素都相等，则返回全 0 向量，避免除零

vmin = min(v);
vmax = max(v);

if abs(vmax - vmin) < eps
    v_n = zeros(size(v));
else
    v_n = (v - vmin) ./ (vmax - vmin);
end
end