function [x, fval, exitflag, output, grad, B] = bfgs(error_func, x0, opt)
% Simple BFGS minimization.
% TEST: David's new BFGS regularization and conditioning implementation for Spinach.

% Ville Bergholm 2016


%% default params

if 1
ls_params = struct(...
    'method', 'bracket-section',...
    'rules', 'Wolfe-strong',...
    'tols', struct(...
        'c1', 1e-2,...  % func
        'c2', 0.9,...   % gradient
        'tau1', 3,...   % bracket expansion
        'tau2', 0.1,...   % left bracket reduction
        'tau3', 0.5...   % right bracket reduction
        ));
elseif 1
ls_params = struct(...
    'method', 'backtracking',...
    'rules', 'Armijo',...
    'tols', struct(...
        'c1', 1e-2,...
        'tau', 2/(1+sqrt(5))...  % step reduction factor
        ));
else
ls_params = struct(...
    'method', 'backtracking',...
    'rules', 'Goldstein',...
    'tols', struct(...
        'c1', 0.25,...
        'tau', 2/(1+sqrt(5))...
        ));
end

reg_params = struct(...
    'method', 'RFO',...             % '', 'RFO','TRM','CHOL'
    'cond_method', 'iterative',...  % 'iterative','scaled'
    'n_reg', 2500,...               % max iterates
    'alpha', 1,...
    'delta', 1,...
    'phi', 0.9,...
    'max_cond', 1e5,...
    'track_eig', false);

% use the same field names as fminunc to make comparison easy.
opt_defaults = struct(...
    'MaxIterations', 1000,...
    'MaxFunctionEvaluations', 1000,...
    'FunctionTolerance', NaN,...  % error func target
    'OptimalityTolerance', 1e-6,...
    'StepTolerance', 1e-6,...
    'Display', 'final',...
    'OutputFcn', []);

opt = apply_options(opt_defaults, opt, true);

% initial error and gradient
x = x0(:);
[fval, grad] = error_func(x);

exitflag = [];

output.fval = fval;
output.n_iter = 0;
output.n_evals = 1;
output.stepsize = NaN;
output.grad_norm = NaN;
output.message = [];

% initial inverse Hessian approximation
n_vars = length(x);
B = eye(n_vars);

% iteration loop
while isempty(exitflag)

    %% check, regularize and/or condition the inverse Hessian

    % tests for positive definite and initialise exitflag
    [~, p] = chol(B);
    if ~issymmetric(B), exitflag = -4;
    elseif ~isreal(B), exitflag =-5;
    elseif (~isempty(reg_params.method) && p) ||...
            (~isempty(reg_params.cond_method) && cond(B) > reg_params.max_cond)

        disp('Hessian regularization')
        [B, grad, output] = hessreg(B, grad, reg_params, output);
    elseif reg_params.track_eig
        % explicit tracking of Hessian eigenvalues (for debugging)
        [v,d] = eig(B);
        output.hess_eigvecs = v;
        output.hess_eigs = diag(d);
        output.hess_mineig = min(output.hess_eigs);
        output.hess_cond = max(output.hess_eigs) / output.hess_mineig;
        output.n_cond = 0;
    end

    % final positive definite test
    [~, p] = chol(B);
    if p, exitflag = -3; end

    %% find the search direction

    dir = -B * grad;

    %% call output function

    if ~isempty(opt.OutputFcn) && opt.OutputFcn(x, output, 'iter'), exitflag = -1; end

    if ~isempty(exitflag)
        break
    end

    %% make a line search in the search direction
    % also computes fval and grad for the next point

    % remember current gradient
    old_grad = grad;

    [alpha, fval, grad, exitflag, output] = linesearch(error_func, dir, x, fval, grad, output, ls_params);
    % NOTE that linesearch computes fval and grad at the new point

    output.fval = fval;
    output.grad_norm = norm(grad, Inf);  % first-order optimality measure

    if ~isempty(exitflag)
        break
    end

    %% update current x with search direction and step length

    step = alpha*dir;  % alpha == 1 would be the Newton step
    x = x +step;

    output.stepsize = max(abs(step));  % == norm(step, Inf)

    %% Hessian update

    % gradient difference
    y = grad -old_grad;

    % BFGS update for the inverse Hessian
    temp0 = step' * y;  % scalar
    temp1 = B * y;      % vector
    temp2 = (temp0 +y'*temp1) / temp0^2;  % scalar
    temp3 = temp1 * step';  % rank-1 matrix
    B_new = B +(temp2 * step) * step' -(temp3 +temp3') / temp0;

    % eliminate rounding errors
    B = (B_new +B_new') / 2;


    %% iteration done

    output.n_iter = output.n_iter+1;

    %% check termination conditions (note order, more important at the end)

    if output.n_iter >= opt.MaxIterations || output.n_evals >= opt.MaxFunctionEvaluations,  exitflag = 0; end
    if fval             < opt.FunctionTolerance,   exitflag = 4; end
    % TODO actualerror change, exitflag = 3
    if output.stepsize  < opt.StepTolerance,       exitflag = 2; end
    if output.grad_norm < opt.OptimalityTolerance, exitflag = 1; end
end


% Call output function if specified
if ~isempty(opt.OutputFcn) && opt.OutputFcn(x, output, 'done'), exitflag = -1; end

% exiflag message
switch exitflag
  case  0, temp='Number of iterations/evaluations exceeded MaxIterations/MaxFunctionEvaluations.';
  case  1, temp='Magnitude of gradient smaller than OptimalityTolerance.';
  case  2, temp='Change in x smaller than StepTolerance.';
  case  3, temp='Change in error function less than FunctionTolerance.';
  case  4, temp='Error function value FunctionTolerance reached.';
  case -1, temp='Algorithm was terminated by the output function.';
  case -2, temp='Line search cannot find an acceptable point along the current line.';
  case -3, temp='Hessian matrix is not positive definite.';
  case -4, temp='Hessian matrix is not symmetric.';
  case -5, temp='Hessian matrix is not real.';
  otherwise, error('Undefined exit code.');
end
output.message = temp;
end



function res = issymmetric(A)
% Returns true iff the matrix A is symmetric.
    res = norm(A-A.') < 1e-10;
end
