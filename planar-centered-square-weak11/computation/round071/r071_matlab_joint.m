function r071_matlab_joint(inputDir, augmentedMat, outputPath)
% Independent floating LP solve. Exact certification is performed in Python.
% Reproduce after archive from computation/round071:
% r071_matlab_joint('../inputs/r071','../results/round071-joint/one-column.mat', ...
%     '../results/round071-root/matlab-joint.json')
if nargin < 3
    error('Supply an input directory, augmented MAT path and output JSON path.');
end
cases = {'six-active', 'single-positive', 'one-column'};
inputs = {fullfile(inputDir, 'six-active.mat'), ...
    fullfile(inputDir, 'single-positive.mat'), augmentedMat};
expected = [-1/100, 7/100, 101/10000];
result.schema = 'r071-matlab-joint-numeric-v1';
result.matlab = version;
result.pid = feature('getpid');
result.optimization_toolbox = ver('optim');
result.scope = ['Independent floating solves of three exported finite LPs; ', ...
    'not an exact certificate or a continuous-position result.'];
options = optimoptions('linprog', 'Display', 'none', ...
    'ConstraintTolerance', 1e-9, 'OptimalityTolerance', 1e-9);
result.options = struct('ConstraintTolerance', 1e-9, ...
    'OptimalityTolerance', 1e-9, 'Algorithm', options.Algorithm);
for k = 1:numel(cases)
    data = load(inputs{k});
    data.b = data.b(:);
    data.f = data.f(:);
    data.lb = data.lb(:);
    data.ub = data.ub(:);
    assert(size(data.A, 2) == 145 + (k == 3));
    assert(data.f(end) == -1);
    started = tic;
    [x, fval, exitflag, output, lambda] = linprog(data.f, data.A, ...
        data.b, [], [], data.lb, data.ub, options);
    item.name = cases{k};
    item.elapsed_seconds = toc(started);
    item.exitflag = exitflag;
    item.output = output;
    item.rows = size(data.A, 1);
    item.columns = size(data.A, 2);
    item.nonzeros = nnz(data.A);
    assert(exitflag > 0, 'linprog did not report success.');
    item.objective = fval;
    item.eta = x(end);
    item.expected_eta = expected(k);
    item.max_constraint_violation = max([0; data.A*x-data.b]);
    item.min_outside_mass = min(x(1:end-1));
    item.primal = x;
    nz = find(lambda.ineqlin ~= 0);
    item.nonzero_inequality_dual = [nz, lambda.ineqlin(nz)];
    assert(abs(item.eta-expected(k)) < 1e-8);
    assert(item.max_constraint_violation < 1e-7);
    assert(item.min_outside_mass > -1e-8);
    result.cases(k) = item;
end
result.numeric_checks_pass = true;
folder = fileparts(outputPath);
if ~exist(folder, 'dir')
    mkdir(folder);
end
fid = fopen(outputPath, 'w');
assert(fid >= 0, 'Cannot open result output.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(result, PrettyPrint=true));
fprintf('R71 MATLAB independent numeric checks PASS; PID %d\n', result.pid);
for k = 1:numel(cases)
    fprintf('%s: eta=%.16g, max violation=%.3g\n', ...
        result.cases(k).name, result.cases(k).eta, ...
        result.cases(k).max_constraint_violation);
end
end
