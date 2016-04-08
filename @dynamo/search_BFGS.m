function [exitflag, output] = search_BFGS(self, obj_func, matlab_options)
% BFGS optimization.

% define the optimization problem
problem.objective = obj_func;
problem.x0 = self.seq.get(self.opt.control_mask);
problem.solver = 'fminunc';

% default options for fminunc

% TODO with newer MATLAB versions we would do it like this:
%problem.options = optimoptions('fminunc',...
%    'Algorithm',    'quasi-newton',...
% for now, use the old optimset()
problem.options = optimset(...
    'DerivativeCheck', 'off',...
    'Display',         'final',...
    'GradObj',         'on',...    % use user-supplied gradient
    'LargeScale',      'off', ...  % force quasi-newton algorithm (BFGS)
    'MaxIter',         1e4,...
    'OutputFcn', @(x, optimValues, state) monitor_func(self, x, optimValues, state),...
    'TolFun',          1e-8,...
    'TolX',            1e-8);
% additional user-defined options
problem.options = optimset(problem.options, matlab_options);
% save a copy
self.opt.matlab_options = problem.options;

fprintf('\nOptimizing algorithm: BFGS. Running...\n\n'); drawnow;

% try to minimise objective function to zero
[x, cost, exitflag, output] = fminunc(problem);

% minimizer may be different than the last point evaluated
self.update_controls(x, self.opt.control_mask);
end
