function result = propose_mass_joint(inputPath, outputPath)
% Heuristic LPs only; every output requires exact feasibility and halo checks.
data = jsondecode(fileread(inputPath));
W = data.W;
G = data.G;
n = size(W, 2);
assert(size(W, 1) == numel(data.mass));
assert(size(G, 2) == n);
identity = eye(n);
options = optimoptions('linprog', 'Display', 'none');
proposals = struct([]);
started = tic;
for gi = 1:size(G, 1)
    g = G(gi, :)';
    normalization = max(abs(g));
    assert(normalization > 0);
    f = [-g / normalization; zeros(n, 1)];
    for ti = 1:numel(data.trust)
        sigma = data.trust(ti);
        for bi = 1:numel(data.budgets)
            budget = data.budgets(bi);
            A = [identity, -identity; -identity, -identity; ...
                 zeros(1, n), ones(1, n); ...
                 sigma * W, zeros(size(W)); ...
                 -sigma * W, zeros(size(W))];
            b = [zeros(2*n, 1); budget; ...
                 data.upper(:) - data.mass(:); data.mass(:) - data.lower(:)];
            lb = [-ones(n, 1); zeros(n, 1)];
            ub = ones(2*n, 1);
            [x, fval, flag, details] = linprog(f, A, b, [], [], lb, ub, options);
            assert(flag > 0, 'LP did not finish successfully');
            p.gradient_index = gi;
            p.gradient_scale = data.gradient_scales{gi};
            p.trust_index = ti;
            p.trust_exact = data.trust_exact{ti};
            p.budget = budget;
            p.u = x(1:n)';
            p.exitflag = flag;
            p.iterations = details.iterations;
            p.predicted_gain = sigma * g' * x(1:n);
            p.scaled_objective = fval;
            p.max_constraint_violation = max([0; A*x-b; lb-x; x-ub]);
            if isempty(proposals)
                proposals = p;
            else
                proposals(end+1) = p; %#ok<AGROW>
            end
        end
    end
end
result.schema = 'matlab-mass-lp-proposals-v1';
result.scope = 'Heuristic proposals only; no derivative or bound certificate';
result.version = version;
result.pid = feature('getpid');
result.exact_input_sha256 = data.exact_input_sha256;
result.wall_seconds = toc(started);
result.lp_count = numel(proposals);
result.proposals = proposals;
fid = fopen(outputPath, 'w');
assert(fid >= 0);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(result, PrettyPrint=true));
end
