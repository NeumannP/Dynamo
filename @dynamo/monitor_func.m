function stop = monitor_func(self, x, optimValues, state)
% Executed once every iteration during optimization, decides if we should stop here.

stop = false;
wt = (now() -self.opt.wall_start) * 24*60*60; % elapsed time in seconds
ct = cputime() -self.opt.cpu_start;


%% Check for user interrupt

drawnow(); % flush drawing events and callbacks (incl. user interrupts)

% check for stop signal from the UI
% TODO this is where we would check for Ctrl-C if MATLAB supported it...
if self.config.stop
    self.opt.term_reason = 'User interrupt';
    stop = true;
end


%% Plot the sequence every now and then

% Refresh the UI figure.
temp = self.opt.options.plot_interval;
if ~isempty(self.config.ui_fig) && temp && mod(self.opt.n_iter, temp) == 0
    self.ui_refresh(false, optimValues.fval);
end


%% Check termination conditions

% TODO some of these are already present in optimValues...
self.opt.n_iter = self.opt.n_iter +1;

if self.opt.n_eval >= self.opt.options.max_evals
    self.opt.term_reason = 'Evaluation limit reached';
    stop = true;
end

if wt >= self.opt.options.max_walltime
    self.opt.term_reason = 'Wall time limit reached';
    stop = true;
end

if ct >= self.opt.options.max_cputime
    self.opt.term_reason = 'CPU time limit reached';
    stop = true;
end

if self.opt.last_grad_norm <= self.opt.options.min_gradient_norm
    self.opt.term_reason = 'Minimum gradient norm reached';
    stop = true;
end

% have we reached our goal?
if optimValues.fval <= self.opt.options.error_goal
    self.opt.term_reason = 'Error goal achieved';
    stop = true;
end


%% Stats collector part

self.opt.error(end+1) = optimValues.fval;
self.opt.wall_time(end+1) = wt;
self.opt.cpu_time(end+1)  = ct;
self.opt.control_integral(end+1,:) = self.seq.integral();
end
