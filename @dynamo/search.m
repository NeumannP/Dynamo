function term_reason = search(self, control_mask, varargin)
% Run the optimization.
% If control_mask is empty, try to continue the previous optimization.

if nargin < 2
    % by default update all the controls except taus
    control_mask = self.full_mask(false);
end

%% MATLAB-style options processing.
% Converts a list of fieldname, value pairs in varargin to a struct.
user_options = struct(varargin{:});

if isempty(control_mask)
    % continue previous optimization (re-use most of the opt struct)
    control_mask = self.opt.control_mask;
    % self.opt.control_mask, self.opt.options and self.opt.matlab_options are kept as is
else
    % default termination conditions and other options
    self.opt.options = struct(...
        'error_goal',        0.5 * (1e-4)^2 / self.system.norm2,...
        'max_evals',         1e6,...
        'max_walltime',      1800,...
        'max_cputime',       1e6,...
        'min_gradient_norm', 1e-20,...
        'plot_interval',     1);   % how often should we plot intermediate results?
    self.opt.matlab_options = struct();
end
% initialization of the optimization data structures
self.init_opt();
self.opt.control_mask = control_mask;

% modify self.opt.options with user_options (if any)
[self.opt.options, unused] = apply_options(self.opt.options, user_options, true);
% the rest are dumped into matlab_options
[self.opt.matlab_options] = apply_options(self.opt.matlab_options, unused, false);


%% run the optimizer

fprintf('Optimization space dimension: %d\n', sum(sum(self.opt.control_mask)));

% define the optimization problem
obj_func = @(x) goal_and_gradient_function_wrapper(self, x);

% run BFGS optimization
[self.opt.matlab_exitflag, self.opt.matlab_output] = self.search_BFGS(obj_func, self.opt.matlab_options);


%% make a copy of the current run's statistics

self.stats{end+1} = self.opt;


%% termination reason

if self.opt.matlab_exitflag == -1
    % term_reason comes from monitor_func
else
    % terminated by optimizer
    self.opt.term_reason = self.opt.matlab_output.message;
end
term_reason = self.opt.term_reason;
end


function [err, grad] = goal_and_gradient_function_wrapper(self, x)
% x is a vector containing (a subset of) the controls

    self.opt.n_eval = self.opt.n_eval +1;

    self.update_controls(x, self.opt.control_mask);
    [err, grad] = self.compute_error(self.opt.control_mask);
    self.opt.last_grad_norm = sqrt(sum(sum(grad .* grad)));
end
