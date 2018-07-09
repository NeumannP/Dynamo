function term_reason = search_NM(self, control_mask, varargin)
% Nelder-Mead simplex method optimization.


% common initialization of the optimization data structures
user_options = self.init_opt(control_mask, varargin{:});

% define the optimization problem
problem.objective = @(x) goal_function_wrapper(self, x);
problem.x0 = self.seq.get_raw(self.opt.control_mask);
problem.solver = 'fminsearch';

% default options for fminsearch
problem.options = optimset(...
    'Display',         'final',...
    'MaxIter',         1e4,...
    'OutputFcn', @(x, optimValues, state) monitor_func(self, x, optimValues, state),...
    'TolFun',          1e-3,...
    'TolX',            1e-3);
% additional user-defined options
problem.options = optimset(problem.options, user_options);
% save a copy
self.opt.NM_options = problem.options;

fprintf('\nOptimizing algorithm: Nelder-Mead. Running...\n\n'); drawnow;

% try to minimise objective function to zero
[x, cost, exitflag, output] = fminsearch(problem);

self.update_controls(x, self.opt.control_mask); % It may be different than the last point evaluated
term_reason = self.opt.term_reason;
end


function [err] = goal_function_wrapper(self, x)
% x is a vector containing (a subset of) the controls
    
    self.opt.N_eval = self.opt.N_eval + 1;
    self.update_controls(x, self.opt.control_mask);
    [err] = self.compute_error();
end
