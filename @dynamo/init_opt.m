function init_opt(self)
% Initialize the run-specific optimization data and statistics.

%% statistics

self.opt.error = self.compute_error();
self.opt.wall_time = 0;
self.opt.cpu_time  = 0;
self.opt.control_integral  = self.seq.integral();


%% other optimization data

self.opt.control_mask = [];
self.opt.initial_controls = self.seq.get_raw();
self.opt.n_iter = 0;
self.opt.n_eval = 0;
self.opt.wall_start = now();
self.opt.cpu_start = cputime();
self.opt.last_grad_norm = NaN;
self.opt.max_violation = 0;  % track the worst gradient approximation violation
self.opt.term_reason = 'None yet';

self.config.stop = false;  % communication between the UI figure and monitor_func
end
