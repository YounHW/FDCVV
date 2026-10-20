%% MTE_function.m
% MTE = MEPE (CV evidence + variance) + TEAD Taylor nonlinearity enhancement
% Candidate-enumeration version (沿用 TEAD 的 lhs_scaled + knnsearch 框架)
%
% Interface: x_new = MTE_function(obj, A, strategy)
function x_new = MTE_function(obj, A, strategy)

    % =========================================================
    % 0) 读取/设置可选参数（与 TEAD/MEPE 接口保持兼容）
    % =========================================================
    beta = 0.5;      % 增强项权重：RC_MTE = (MEPE_norm) + beta * w * (Res_norm)
    epsilon = 0.01;     % 非聚簇硬约束阈值：若 Dmin(x) < epsilon 则丢弃该候选点

    % ===== 读取消融开关 =====
    use_MEPE = 1;
    use_Taylor = 1;
    use_gate = 1;
    use_epsilon = 1;
    use_adaptive_alpha = 1;
    alpha_fixed = 0.5;

    if isstruct(strategy)
        if isfield(strategy, 'beta');    beta = strategy.beta; end
        if isfield(strategy, 'epsilon'); epsilon = strategy.epsilon; end
    end

    % =========================================================
    % 1) 与 MEPE_function 一致：用 MEPE_q / old_Model 计算 alpha
    %    global 变量由外层循环维护（你需要确保 MEPE_q 在外层正确初始化）
    % =========================================================
    global MEPE_q
    global old_Model

    if MEPE_q == 1
        % 第一次调用：MEPE_function 中 e_true = inf
        old_Model = obj;
        e_true = inf;
    elseif MEPE_q > 1
        % 后续调用：计算“上一轮模型”对最新样本点的误差
        e_true = abs(old_Model.predict(obj.X(end,:)) - obj.Y(end));
        old_Model = obj;
    end

    % =========================================================
    % 2) 复用 MEPE_function 中的：
    %    (1) compute_beta_hat / compute_sigma_sq_hat
    %    (2) sigma_Y_sq_hat_fun(x) = compute_sigma_Y_sq_hat(...)
    % =========================================================
    beta_hat   = compute_beta_hat(obj, obj.R);
    sigma_sq_hat = compute_sigma_sq_hat(obj, obj.R, beta_hat);

    % sigma_Y_sq_hat_fun(x): MEPE 中的方差项（预测误差方差部分）
    sigma_Y_sq_hat_fun = @(x) compute_sigma_Y_sq_hat(obj, sigma_sq_hat, x, obj.theta, obj.R);

    % =========================================================
    % 3) 复用 MEPE_function 中的 LOOCV 快速误差 e_CV_fast
    %    e_CV(i) 对应训练点 obj.X(i,:) 处的 LOOCV evidence
    % =========================================================
    e_CV = e_CV_fast_function(obj);   % obj.m x 1

    % MEPE 的 alpha（与 MEPE_alpha_function 完全一致）
    % alpha = MEPE_alpha_function(obj, e_true, e_CV(end));
    if use_adaptive_alpha == 1
        alpha = MEPE_alpha_function(obj, e_true, e_CV(end));
    else
        alpha = alpha_fixed;
    end

    % =========================================================
    % 4) TEAD 的梯度计算（中心差分），delta_s(i,:) = grad M_hat(x_i)
    %    该梯度用于 Taylor 一阶近似
    % =========================================================
    delta_s = central_difference_gradient(obj); % obj.m x n

    % =========================================================
    % 5) 生成候选点（沿用 TEAD_function）
    % =========================================================
    n = size(obj.X, 2);
    lb = A(1,:);
    ub = A(2,:);

    number_candidate_points = 20000 * n;
    x_cand = lhs_scaled(number_candidate_points, lb, ub); % N x n

    % =========================================================
    % 6) 对每个候选点找到最近训练点：
    %    - Idx(k) 给出最近训练点索引（i*(x)）
    %    - D(k) 给出最近距离（Dmin(x)）
    % =========================================================
    [Idx, D] = knnsearch(obj.X, x_cand); % N x 1 each

    % =========================================================
    % 7) TEAD 的距离门控权重 w(x)
    %    TEAD_function: Lmax = max(max(pdist2(obj.X,obj.X)));
    %    w = 1 - D/Lmax
    % =========================================================
    Lmax = max(max(pdist2(obj.X, obj.X)));
    w = ones(number_candidate_points, 1) - D ./ Lmax;
    if use_gate == 0
        w = ones(number_candidate_points,1);
    end

    % =========================================================
    % 8) 在候选点上计算两部分：
    %    (A) MEPE acquisition：alpha * e_CV(nearest) + (1-alpha) * sigma_Y_sq_hat_fun(x)
    %    (B) TEAD Taylor residual：
    %        g = x_i - x
    %        t(x) = Y_i + grad_i^T * (x_i - x)
    %        s(x) = obj.predict(x)
    %        Res(x) = norm(s(x) - t(x))
    %
    %    然后组合：
    %      J_MTE = MEPE_norm + beta * w * Res_norm
    %
    %    同时加入可选非聚簇硬约束：若 D < epsilon 则丢弃
    % =========================================================
    N = number_candidate_points;

    % 预分配数组（提升速度，避免在循环中动态扩展）
    MEPE_val = -inf(N,1);
    Res = zeros(N,1);

    % （可选）如果 epsilon>0，会丢弃接近已有样本的候选
    if use_epsilon == 1 && epsilon > 0
        validMask = (D >= epsilon);
    else
        validMask = true(N,1);
    end

    for k = 1:N
        if ~validMask(k)
            % 丢弃该候选点：MEPE_val 保持 -inf，则最终不会被选中
            continue;
        end

        % 最近训练点索引 i*(x)
        iStar = Idx(k,1);

        % -------------------------
        % (A) MEPE acquisition
        % -------------------------
        % 注意：你的 MEPE_function 中 fun(x) = -( alpha * e_CV_nearest + (1-alpha)*sigma_Y_sq_hat_fun(x) )
        % 其中 e_CV_nearestPoint_function 直接返回 e_CV(k)（未取 abs/平方，alpha_function 内含 e_true/e_CV 的平方比例）
        % 因此这里严格照抄同样的组合：alpha * e_CV(iStar) + (1-alpha)*sigma_Y_sq_hat_fun(x)
        xk = x_cand(k,:); % 1 x n

%         e_cv_near = e_CV(iStar);                 % LOOCV evidence at nearest training point
%         sigma_part = sigma_Y_sq_hat_fun(xk);   % variance term
%         MEPE_val(k) = alpha * e_cv_near + (1 - alpha) * sigma_part;

        if use_MEPE == 1
            e_cv_near = e_CV(iStar);
            sigma_part = sigma_Y_sq_hat_fun(xk);
            MEPE_val(k) = alpha * e_cv_near + (1 - alpha) * sigma_part;
        else
            MEPE_val(k) = 0;
        end

        % -------------------------
        % (B) TEAD Taylor residual
        % -------------------------
        % TEAD_function 的 Taylor expansion:
        %   g = obj.X(iStar,:) - x_cand(k,:)
        %   t = obj.Y(iStar) + delta_s(iStar,:)*g'
        %   s = obj.predict(x_cand(k,:))
        %   Res = norm(s - t)
%         xi = obj.X(iStar,:);     % 1 x n
%         g = xi - xk;             % 1 x n
% 
%         t = obj.Y(iStar) + delta_s(iStar,:) * g'; % 标量（delta_s * g'）
%         s = obj.predict(xk);                         % 标量预测均值
%         Res(k) = norm(s - t);                       % 与 TEAD_function 一致（虽然是标量，但仍用 norm）

        if use_Taylor == 1
            xi = obj.X(iStar,:);
            g = xi - xk;
            t = obj.Y(iStar) + delta_s(iStar,:) * g';
            s = obj.predict(xk);
            Res(k) = norm(s - t);
        else
            Res(k) = 0;
        end
    end

    % =========================================================
    % 9) 对 MEPE 和 Res 做归一化（用于稳定混合 beta）
    %    这样 beta 的量纲就更好调参，也能避免 MEPE 与 Res 尺度差太大。
    % =========================================================
    % MEPE 归一化：与 TEAD 一样偏向“相对尺度”比较
    if use_MEPE == 1
        validMEPE = isfinite(MEPE_val);
        if any(validMEPE)
            maxMEPE = max(MEPE_val(validMEPE));
            if maxMEPE <= 0 || isnan(maxMEPE)
                maxMEPE = 1;
            end
            MEPE_norm = MEPE_val / maxMEPE;
        else
            MEPE_norm = zeros(N,1);
        end
    else
        MEPE_norm = zeros(N,1);
    end


    % Res 归一化
    if use_Taylor == 1
        maxRes = max(Res);
        if maxRes <= 0 || isnan(maxRes)
            maxRes = 1;
        end
        Res_norm = Res / maxRes;
    else
        Res_norm = zeros(N,1);
    end

    % =========================================================
    % 10) 最终 MTE 打分（MTE核心准则）
    %     J_MTE(x) = MEPE_norm(x) + beta * w(x) * Res_norm(x)
    % =========================================================
    J_MTE = MEPE_norm + beta * (w .* Res_norm);

    % 被丢弃的点（MEPE_val=-inf）会自然不可能取得最大值
    [~, idx_max] = max(J_MTE);

    % =========================================================
    % 11) 返回下一采样点
    % =========================================================
    x_new = x_cand(idx_max,:);

end


% =========================================================
% ========== 局部函数：严格复用你的 MEPE/TEAD 计算方式 ==========
% =========================================================

function e_CV_fast = e_CV_fast_function(obj)
    % MEPE_function.m 里的 e_CV_fast_function（原样一致）
    beta_hat = compute_beta_hat(obj, obj.R);
    H = obj.F *((obj.F' * obj.F) \ obj.F');
    d = obj.Y - obj.F * beta_hat;
    inv_R = inv(obj.R);

    e_CV_fast = zeros(obj.m,1);
    for i=1:obj.m
        e_CV_fast(i) = ((inv_R(i,:) * ( d + H(:,i) * (d(i)/(1-H(i,i)))) ) ...
                        / inv_R(i,i))^2;
    end
end

function alpha = MEPE_alpha_function(obj, e_true, e_CV)
    % MEPE_function.m 里的 MEPE_alpha_function（原样一致）
    global MEPE_q
    if MEPE_q == 1
        alpha = 0.5;
    elseif MEPE_q > 1
        val = 0.5 * ((e_true^2) / e_CV);
        alpha = 0.99 * min(val, 1);
    end
    MEPE_q = MEPE_q + 1;
end

function delta_s = central_difference_gradient(obj)
    % TEAD_function.m 里的 central_difference_gradient（原样一致）
    h = 1e-4;
    n = size(obj.X,2);
    delta_s = zeros(obj.m,n);
    for i=1:obj.m
        X_ip1 = obj.X(i,:) + h;
        X_im1 = obj.X(i,:) - h;

        Y_ip1 = obj.predict(X_ip1);
        Y_im1 = obj.predict(X_im1);

        delta_s(i,:) = (Y_ip1 - Y_im1) / (2*h);
    end
end