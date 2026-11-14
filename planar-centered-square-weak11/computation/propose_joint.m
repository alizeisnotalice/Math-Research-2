function result = propose_joint(inputPath, outputPath)
% Numerical proposal only. Exact halo verification must certify each input.
data = jsondecode(fileread(inputPath));
G = data.gradients;
kinds = string(data.kinds(:));
n = numel(kinds);
assert(size(G, 2) == n);
groups = ["all", "mass", "coordinate", "geometry"];
budgets = data.budgets(:)';
A = [speye(n), -speye(n); -speye(n), -speye(n); sparse(1,n), sparse(ones(1,n))];
options = optimoptions('linprog', 'Display', 'none');
records = struct('gradient_row', {}, 'group', {}, 'budget', {}, 'u', {}, ...
    'exitflag', {}, 'predicted_gain', {}, 'iterations', {});
started = tic;
for gi = 1:size(G,1)
    gradient = G(gi,:)';
    scale = max(abs(gradient));
    if scale == 0
        scale = 1;
    end
    for group = groups
        active = true(n,1);
        if group == "mass"
            active = kinds == "mass";
        elseif group == "coordinate"
            active = kinds == "coordinate";
        elseif group == "geometry"
            active = kinds ~= "mass";
        end
        lower = [-double(active); zeros(n,1)];
        upper = [double(active); double(active)];
        for budget = budgets
            rhs = [zeros(2*n,1); budget];
            [solution, ~, flag, output] = linprog([-gradient/scale; zeros(n,1)], ...
                A, rhs, [], [], lower, upper, options);
            assert(flag > 0, 'Proposal LP failed: exitflag %d', flag);
            rec.gradient_row = gi;
            rec.group = char(group);
            rec.budget = budget;
            rec.u = solution(1:n)';
            rec.exitflag = flag;
            rec.predicted_gain = gradient' * solution(1:n);
            rec.iterations = output.iterations;
            records(end+1) = rec; %#ok<AGROW>
        end
    end
end
result = struct('scope', 'finite difference LP proposals; no derivative or global optimality certificate', ...
    'matlab_version', version, 'parameters', n, 'lp_count', numel(records), ...
    'elapsed_seconds', toc(started), 'proposals', records);
fid = fopen(outputPath, 'w');
assert(fid >= 0, 'Cannot open proposal output');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(result));
fprintf('proposal_count=%d elapsed_seconds=%.6f\n', numel(records), result.elapsed_seconds);
end
