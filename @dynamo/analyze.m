function analyze(self)
% Analyzes the results of an optimization run.


fprintf('Final normalized error: %g\n    Wall time: %g s\n    CPU  time: %g s\nTermination reason: %s\n\n\n', ...
	self.opt.error(end), self.opt.wall_time(end), self.opt.cpu_time(end), self.opt.term_reason);

fprintf('Number of gradient evaluations: %d\n', self.opt.n_eval);
fprintf('Final sequence duration: %g\n', sum(self.seq.tau));


% rough error estimates
%d = real(eig(self.system.A));
%n = length(d);
%T = sum(self.seq.tau);
%e_max = 1-exp(-T*sum(d)/n)
%e_min = 1-sum(exp(-T*d))/n

% plot the final sequence and some analytics
figure()
ax = subplot(3, 1, 1);
self.plot_seq(ax);


%% plot the error

ax = subplot(3, 1, 2);
set_plotstyle(ax);
self.plot_stats('error semilog', ax);
title(ax, 'Optimization error')


%% plot control integrals

ax = subplot(3, 1, 3);
set_plotstyle(ax);
%set(ax, 'LineStyleOrder','--')
self.plot_stats('control_integral', ax);
title(ax, 'Control integral')
end
