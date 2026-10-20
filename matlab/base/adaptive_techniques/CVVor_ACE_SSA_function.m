function [x_new, details] = CVVor_ACE_SSA_function(obj, A, mode)
% CVVor -> ACE/SSA (two-stage; Stage2 restricted to Cand=C{k,2})
%
% mode: 'ACE' or 'SSA'
% opts: optional struct
%   opts.tolDup    (default 1e-12) remove duplicate candidates
%   opts.eps0      (default 0.01) SSA min-distance feasibility threshold
%   opts.SSA_epsMin(default 1e-6) minimal eps in relaxation
%   opts.maxRelax  (default 20)    relaxation steps
%
% Output:
%   x_new selected point in Cand.

if nargin < 4, opts = struct(); end
if nargin < 3 || isempty(mode), error('mode required'); end
mode = upper(mode);

if ~any(strcmp(mode, {'ACE','SSA'}))
    error("mode must be 'ACE' or 'SSA'.");
end

if ~isfield(opts,'tolDup'),     opts.tolDup = 1e-12; end
if ~isfield(opts,'eps0'),       opts.eps0 = 0.01; end
if ~isfield(opts,'SSA_epsMin'), opts.SSA_epsMin = 1e-6; end
if ~isfield(opts,'maxRelax'),   opts.maxRelax = 20; end



lb = A(1,:); ub = A(2,:);

% ===================== Stage 1: CVVor region =====================
C = randomVoronoi(obj.X, lb, ub);              % your function
eLOO_perPoint = eL00_function(obj);            % same as your function

[~, k_star] = max(eLOO_perPoint);
Cand = C{k_star,2};                            % candidates in sensitive Voronoi cell

% remove duplicates to existing samples
if ~isempty(Cand)
    keep = true(size(Cand,1),1);
    for i = 1:size(Cand,1)
        if min(vecnorm(obj.X - Cand(i,:), 2, 2)) <= opts.tolDup
            keep(i) = false;
        end
    end
    Cand = Cand(keep,:);
end
if isempty(Cand)
    error('Cand is empty after duplicate removal. Check randomVoronoi settings.');
end

% ===================== Stage 2: Point selection =====================
switch mode
    case 'ACE'
        [alpha, d] = ACE_alpha_local(obj);

        % ACE feasible: min_i ||x - Xi|| >= d
        feasibleMask = true(size(Cand,1),1);
        for j = 1:size(Cand,1)
            feasibleMask(j) = (min(vecnorm(obj.X - Cand(j,:), 2, 2)) >= d - 1e-12);
        end
        Cand_feas = Cand(feasibleMask,:);
        if isempty(Cand_feas)
            Cand_feas = Cand; % fallback if too strict
        end

        % ACE objective: E(x) = sum_i eLOO(i) * exp(-alpha*||Xi-x||)
        % with ACE_optimize convention using E_min = -E, minimize E_min
        % => equivalently maximize E.
        E = zeros(size(Cand_feas,1),1);
        for j = 1:size(Cand_feas,1)
            x = Cand_feas(j,:);
            s = 0;
            for i = 1:obj.m
                s = s + eLOO_perPoint(i) * exp(-alpha * norm(obj.X(i,:) - x));
            end
            E(j) = s;
        end

        [~, idxBest] = max(E);
        x_new = Cand_feas(idxBest,:);

        details = struct();
        details.mode = mode;
        details.k_star = k_star;
        details.Cand_size = size(Cand,1);
        details.Cand_feas_size = size(Cand_feas,1);
        details.alpha = alpha;
        details.d = d;
        details.bestACE_E = E(idxBest);

    case 'SSA'
        % SSA discrete approximation:
        % for feasible candidates x, score(x)=Delta*CMD where
        %   Delta = (mu(x) - mu^{(-p)}(x))^2 , CMD(x)=sum_i ||x-Xi||^2
        epsVal = opts.eps0;

        % compute p order once using CMD on existing points (same spirit as your SSA)
        CMD_existing = zeros(obj.m,1);
        for i=1:obj.m
            val = 0;
            for j=1:obj.m
                val = val + norm(obj.X(i,:) - obj.X(j,:))^2;
            end
            CMD_existing(i) = val;
        end
        [~, p_order] = sort(CMD_existing, 'descend');

        bestScore = -inf;
        x_new = Cand(1,:);
        found = false;

        for relaxIter = 1:opts.maxRelax
            % feasibility: min distance to existing points > epsVal
            feasMask = true(size(Cand,1),1);
            for j=1:size(Cand,1)
                feasMask(j) = (min(vecnorm(obj.X - Cand(j,:), 2, 2)) > epsVal);
            end
            Cand_feas = Cand(feasMask,:);
            if isempty(Cand_feas)
                epsVal = epsVal/2;
                if epsVal < opts.SSA_epsMin
                    break;
                end
                continue;
            end

            % CMD(x) for all feasible candidates
            nF = size(Cand_feas,1);
            CMDx = zeros(nF,1);
            for j = 1:nF
                x = Cand_feas(j,:);
                val = 0;
                for ii=1:obj.m
                    val = val + norm(x - obj.X(ii,:))^2;
                end
                CMDx(j) = val;
            end

            % mu(x) using current model (IMPORTANT: use obj.predict method)
            [mu_hat, ~] = obj.predict(Cand_feas);   % returns vector (nF x 1)

            % iterate pidx in SSA order; build OKp and compute mu^{(-p)}(x)
            for t = 1:obj.m
                pidx = p_order(t);

                Xp = obj.X; Xp(pidx,:) = [];
                Yp = obj.Y; Yp(pidx) = [];

                OKp = OK_model(obj.auto_correlation_function, Xp, Yp, obj.theta_opti_technique);

                [mu_hatp, ~] = OKp.predict(Cand_feas);  % (nF x 1)

                Delta = (mu_hat - mu_hatp).^2;          % (nF x 1)
                score = Delta .* CMDx;                  % (nF x 1)

                [sBestThis, jBestThis] = max(score);
                if sBestThis > bestScore
                    bestScore = sBestThis;
                    x_new = Cand_feas(jBestThis,:);
                    found = true;
                end
            end

            if found
                break;
            end

            epsVal = epsVal/2;
            if epsVal < opts.SSA_epsMin
                break;
            end
        end

        details = struct();
        details.mode = mode;
        details.k_star = k_star;
        details.Cand_size = size(Cand,1);
        details.eps_used = epsVal;
        details.bestScore = bestScore;

end

end


% ===================== your CVVor LOO error =====================
function eLOO_perPoint = eL00_function(obj)
eLOO_perPoint = zeros(obj.m,1);
for i=1:obj.m
    Xp = obj.X; Xp(i,:) = [];
    Yp = obj.Y; Yp(i) = [];

    M_hatp = OK_model(obj.auto_correlation_function, Xp, Yp, obj.theta_opti_technique);

    y = obj.Y(i);
    % IMPORTANT: call class method predict
    [mu_hatp, ~] = predict(M_hatp, obj.X(i,:));
    eLOO_perPoint(i) = norm(y - mu_hatp);
end

end


% ===================== your ACE alpha/d =====================
function [alpha, d] = ACE_alpha_local(obj)
min_distances = zeros(obj.m,1);
for i=1:obj.m
    distance_min = zeros(obj.m-1,1);
    iter = 1;
    for j=1:obj.m
        if i ~= j
            distance_min(iter) = norm(obj.X(i,:) - obj.X(j,:));
            iter = iter + 1;
        end
    end
    min_distances(i) = min(distance_min);
end

max_distance = max(min_distances);
DOI_distance = 0.5 * max_distance;
worst_case_DOI = 1e-5;
alpha = -log(worst_case_DOI) / DOI_distance;

d = 0.5 * mean(min_distances);
end