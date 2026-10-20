function x_new = fCVVor_function(obj,A)
% Cross-Validation-Voronoi (CVVOR) - Fast LOOCV version (方案A)
% Xu, Shengli, et al. (2014)
%
% inputs:
% obj - OK kriging object
% A   - [lb; ub] parametric space bounds, size 2xd
%
% outputs:
% x_new - new found sample point



lb = A(1,:);
ub = A(2,:);

% 1) Build Voronoi cells (discretized)
C = randomVoronoi(obj.X,lb,ub);

% 2) Fast LOOCV error per training point
%    (replace real LOOCV loop)
eLOO_perPoint = eLOO_fast_function(obj);

% 3) choose most sensitive cell (the cell whose representative point is from max LOOCV error)
[~,k] = max(eLOO_perPoint);

% representative point for cell k
P_sensitive = C{k,1};

% 4) Exploration inside that cell: choose farthest candidate from representative point
maxDistance = -inf;
x_new = C{k,2}(1,:); % init (in case all distances fail)
for i=2:size(C{k,2},1)
    distance = norm(C{k,2}(i,:) - P_sensitive);
    if distance > maxDistance
        x_new = C{k,2}(i,:);
        maxDistance = distance;
    end
end

end

% =====================================================================
% Fast LOOCV function
% =====================================================================
function eLOO_perPoint = eLOO_fast_function(obj)
% Returns a vector of size obj.m x 1
% where larger means "more error" (used by argmax)

% ---- Option 1 (recommended): call your existing fast LOOCV function ----
% Change the function name below to your actual implementation.
% Example: e_CV_fast_function(obj) should return obj.m x 1.
try
    eLOO_perPoint = e_CV_fast_function(obj);
catch ME
    error(['Cannot find fast LOOCV function "e_CV_fast_function(obj)". ' ...
           'Please replace it with your real fast LOOCV function name. Error: %s'], ME.message);
end

% ---- Optional: if your fast output is squared error, you may want absolute error. ----
% Your original CVVor uses: norm(y - mu_hatp)  (absolute error for scalar y).
% If e_CV_fast_function returns squared errors, argmax still works, but if you want
% absolute error scaling, uncomment the next line.
% eLOO_perPoint = sqrt(abs(eLOO_perPoint));

% ---- Shape check ----
eLOO_perPoint = eLOO_perPoint(:);
if length(eLOO_perPoint) ~= obj.m
    error('Fast LOOCV output length (%d) does not match obj.m (%d).', length(eLOO_perPoint), obj.m);
end
end
