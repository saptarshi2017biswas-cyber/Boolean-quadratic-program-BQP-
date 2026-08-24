%% Question-2

function Q2()
disp('>>> ENTERED Q2')

%% ================= PARAMETERS =================
n = 100;
rng(1);

%% ================= GENERATE C =================
A = randn(n,n);
C = A * A';   % symmetric PSD matrix

%% ================= SDP RELAXATION =================
cvx_begin sdp quiet
    variable X(n,n) symmetric
    maximize( trace(C*X) )
    subject to
        diag(X) == 1;
        X >= 0;   % PSD constraint
cvx_end

fprintf('SDP optimal value: %.4f\n', cvx_optval);

%% ================= RANDOMIZATION =================
% Cholesky / Eigen decomposition
[V,D] = eig(X);
D = max(D,0);   % numerical safety
X_sqrt = V * sqrt(D);

num_trials = 100;
best_obj = -inf;
best_x = [];

for k = 1:num_trials
    
    % Gaussian random vector
    r = randn(n,1);
    
    % Generate candidate solution
    x_candidate = sign(X_sqrt * r);
    
    % Replace zeros (rare case)
    x_candidate(x_candidate == 0) = 1;
    
    % Evaluate objective
    obj = x_candidate' * C * x_candidate;
    
    if obj > best_obj
        best_obj = obj;
        best_x = x_candidate;
    end
end

fprintf('Best randomized objective: %.4f\n', best_obj);

%% ================= RANK-1 APPROX (OPTIONAL) =================
% Using top eigenvector
[vecs, vals] = eig(X);
[~, idx] = max(diag(vals));
v = vecs(:, idx);

x_rank1 = sign(v);
x_rank1(x_rank1 == 0) = 1;

obj_rank1 = x_rank1' * C * x_rank1;

fprintf('Rank-1 approximation objective: %.4f\n', obj_rank1);

%% ================= COMPARISON =================
fprintf('\n===== COMPARISON =====\n');
fprintf('SDP upper bound        : %.4f\n', cvx_optval);
fprintf('Gaussian randomization : %.4f\n', best_obj);
fprintf('Rank-1 approx          : %.4f\n', obj_rank1);