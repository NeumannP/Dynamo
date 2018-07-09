function [exitflag, output] = search_BFGS(self, obj_func, matlab_options)
% BFGS optimization.


% default options
opt = struct(...
    'MaxIterations',            1e6,...
    'MaxFunctionEvaluations',   1e6,...
    'FunctionTolerance',        NaN,...
    'OptimalityTolerance',     1e-8,...
    'StepTolerance',           1e-8,...
    'Display',              'final');

% apply user-defined options
opt = apply_options(opt, matlab_options, true);
% save a copy
self.opt.matlab_options = opt;

% HACK we do not wish to save the monitor function handle with matlab_options
% because it does not survive saving/loading.
opt.OutputFcn = @(x, optimValues, state) monitor_func(self, x, optimValues, state);

% initial values
x0 = self.seq.get_raw(self.opt.control_mask);

if 0
    % TEST: New BFGS implementation
    fprintf('\nOptimizing algorithm: BFGS. Running...\n\n'); drawnow;
    [x, fval, exitflag, output] = bfgs(obj_func, x0, opt);
else
    % Old BFGS implementation using fminunc
    opt = map_matlab_options(opt, true);

    % additional default options for fminunc
    % TODO with newer MATLAB versions we would do it like this:
    %problem.options = optimoptions('fminunc', 'Algorithm','quasi-newton',...
    opt.DerivativeCheck = 'off';  % no gradient test
    opt.GradObj         =  'on';  % use user-supplied gradient
    opt.LargeScale      = 'off';  % force quasi-newton algorithm (BFGS)

    % define the optimization problem
    problem.solver = 'fminunc';
    problem.objective = obj_func;
    problem.x0 = x0;
    problem.options = opt;

    fprintf('\nOptimizing algorithm: BFGS (fminunc). Running...\n\n'); drawnow;
    [x, fval, exitflag, output] = fminunc(problem);
end

% minimizer may be different than the last point evaluated
self.update_controls(x, self.opt.control_mask);
end


function st = map_matlab_options(st, invert)
% Maps one set of option names to another.
% This is needed for backwards compatibility (for now at least)
% because Matlab now uses optimoptions() as the recommended default
% instead of optimset(), and some of the option names have changed.
% Internally we prefer the newer, more complete options.

% map: {optimset, optimoptions; ...}
    map = {'MaxFunEvals', 'MaxFunctionEvaluations';...
           'MaxIter',     'MaxIterations';...
           'TolFun',      'FunctionTolerance';...
           'TolFun',      'OptimalityTolerance';...
           'TolX',        'StepTolerance';...
          };
    % NOTE: TolFun is mapped to two different fields.

    % inverse mapping?
    if nargin > 1 && invert
        map = fliplr(map);
    end

    % map the fields
    for k=1:size(map, 1)
        temp = map{k, 1}; % field name to be mapped
        if isfield(st, temp)
            st.(map{k, 2}) = st.(temp);
        end
    end
    st = rmfield(st, map(:, 1));  % remove the old fields
end
