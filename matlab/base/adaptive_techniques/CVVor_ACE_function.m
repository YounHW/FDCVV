function x_new = CVVor_ACE_function(obj,A)
% CVVOR-guided ACE: 
% 1) pick the Voronoi cell of the max LOO error point
% 2) select x_new by ACE criterion among discrete candidates in that cell


lb = A(1,:);
ub = A(2,:);

% ----- CVVOR: build Voronoi candidate sets and find most sensitive cell -----
C = randomVoronoi(obj.X,lb,ub);          % C{k,1}=X(k,:), C{k,2}=candidates
eLOO_perPoint = eL00_function(obj);     % LOO error for each existing sample

[~,k] = max(eLOO_perPoint);
candidates = C{k,2};

% ----- ACE: parameters for the ACE criterion and distance constraint -----
[alpha, d] = ACE_alpha(obj);

% ----- Discrete ACE search in the selected Voronoi cell -----
bestScore = -inf;
x_new = candidates(1,:);

for jj = 1:size(candidates,1)
    x = candidates(jj,:);

    % skip if x coincides (or too close) to an existing sample (ACE distance constraint)
    if any(vecnorm(obj.X - x, 2, 2) < d)
        continue
    end

    % ACE score: E(x) = sum e_i * exp(-alpha*||X_i-x||)
    dist = vecnorm(obj.X - x, 2, 2); % m×1
    score = sum(eLOO_perPoint .* exp(-alpha * dist));

    if score > bestScore
        bestScore = score;
        x_new = x;
    end
end

end

% ---------------- LOO error (same as your CVVOR/ACE code) ----------------
function eLOO_perPoint = eL00_function(obj)
eLOO_perPoint = zeros(obj.m,1);
for i=1:obj.m
    Xp = obj.X; Xp(i,:) = [];
    Yp = obj.Y; Yp(i) = [];

    M_hatp = OK_model(obj.auto_correlation_function,Xp,Yp,obj.theta_opti_technique);

    y = obj.Y(i);
    [mu_hatp,~] = predict(M_hatp,obj.X(i,:));

    eLOO_perPoint(i) = norm(y - mu_hatp);
end
end

% ---------------- Reuse your ACE_alpha ----------------
function [alpha, d] = ACE_alpha(obj)
min_distances= zeros(obj.m,1);
for i=1:obj.m
    distance_min = zeros(obj.m-1,1);
    iter = 1;
    for j=1:obj.m
        if ~(i==j)
            distance_min(iter) = norm(obj.X(i,:) - obj.X(j,:));
            iter = iter +1;
        end
    end
    min_distances(i) = min(distance_min);
end

max_distance = max(min_distances);
DOI_distance = 0.5 * max_distance;
worst_case_DOI = 10^(-5);
alpha = -log(worst_case_DOI) / DOI_distance;

d= 0.5 *mean(min_distances);
end