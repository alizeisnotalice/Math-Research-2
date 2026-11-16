function r072_matlab_cells(inputMat, outputPath)
% Numerical re-solve of the final finite cut model, not a geometry certificate.
data = load(inputMat);
f = data.f(:); b = data.b(:); beq = data.beq(:);
lb = data.lb(:); ub = data.ub(:); intcon = double(data.intcon(:));
assert(size(data.A, 2) == 710 && size(data.A, 1) == 8940);
assert(numel(intcon) == 565);
options = optimoptions('intlinprog', 'Display', 'none', ...
    'RelativeGapTolerance', 1e-9, 'AbsoluteGapTolerance', 0);
started = tic;
[x, fval, exitflag, output] = intlinprog(f, intcon, data.A, b, ...
    data.Aeq, beq, lb, ub, options);
result.schema = 'r072-matlab-cell-mip-v1';
result.matlab = version;
result.optimization_toolbox = ver('optim');
result.pid = feature('getpid');
result.elapsed_seconds = toc(started);
result.exitflag = exitflag;
result.output = output;
result.scope = ['Numerical re-solve of the exported final finite cut model. ', ...
    'Unlisted residual geometric constraints are not checked here; ', ...
    'this MATLAB vector is not an additional deletion-area certificate.'];
result.rows = size(data.A, 1); result.columns = size(data.A, 2);
result.integer_variables = numel(intcon);
result.options = struct('RelativeGapTolerance', 1e-9, 'AbsoluteGapTolerance', 0);
assert(~isempty(x), 'No MIP incumbent returned.');
result.objective = fval;
result.producer_objective = data.scipy_objective;
result.objective_difference = fval-data.scipy_objective;
result.cell_objective_divided_by_1e6 = -fval/1e6;
result.max_constraint_violation = max([0; data.A*x-b]);
result.max_bound_violation = max([0; lb-x; x-ub]);
result.max_integer_error = max(abs(x(intcon)-round(x(intcon))));
result.solution = x;
assert(exitflag > 0);
assert(result.max_constraint_violation < 1e-7);
assert(result.max_bound_violation < 1e-7);
assert(result.max_integer_error < 1e-7);
assert(abs(result.objective_difference) < 1e-3);
result.numerical_comparison_pass = true;
fid = fopen(outputPath, 'w'); assert(fid >= 0);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(result, PrettyPrint=true));
fprintf('R72 MATLAB MIP numeric PASS; objective %.12g, PID %d\n', fval, result.pid);
fprintf('Residual %.3g; integer error %.3g; time %.3fs\n', ...
    result.max_constraint_violation, result.max_integer_error, result.elapsed_seconds);
end
