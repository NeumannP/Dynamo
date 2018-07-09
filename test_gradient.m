function [d, direction] = test_gradient(test, d, direction)
% Tests the error functions and their gradients.
%
%  d is a Dynamo instance containing the optimization problem used.
%  If no d is given, uses one of the test suite problems.
%
%  test == 'time':  Measures the walltime required to compute the error function and gradient.
%
%  test == 'acc':  Checks if the computed gradient of the error function is accurate.
%
%  Given an error fuction f(\vec{x}), an accurate gradient
%  evaluated at the point \vec{x0} yields a linearization
%
%    g(\vec{x0} +s \vec{d}) := f(\vec{x0}) +s \vec{d} \dot \vec{grad}.
%
%  \vec{d} is an arbitrary unit vector giving a direction in the parameter space.
%
%  The error |(f-g)(\vec{x0} +s \vec{d})| should scale as O(s^2) (Taylor series).
%  If there is a small error in the gradient, linear scaling
%  will overtake the quadratic one when |s| is small enough:
%
%    |s_cutoff| = 2 |\vec{d} \dot \vec{grad_error}| / |\vec{d}^T Hessian_f(x_0) \vec{d}|
%
%  1st order approximations turn to O(s) scaling almost immediately
%  For basic finite difference methods, \vec{grad_error} should be proportional to epsilon.

% Ville Bergholm 2011-2016


%randseed(seed);

%% set up a system, random controls

if nargin < 2
    switch test
      case 'acc'
        d = test_suite(21);
        d.set_controls([], d.seq.tau+0.01*randn(size(d.seq.tau)));
        
      case 'time'
        d = test_rand_problem('closed gate', 16, 4);
    end
end
   

%% choose an error function and a compatible gradient

ff = 'full'
gg = 'fd'

d.config.epsilon  = 2e-8;
d.config.UL_mixed = false;

ttt = ['error\_', ff, ', gradient\_', gg];

switch gg
  case 'fd'
    ttt = sprintf('%s, epsilon = %g', ttt, d.config.epsilon);
end


d.config.dP = gg;
switch ff
  case 'abs'
    d.config.error_func = @error_abs;
    d.config.nonprojective_error = false;

  case 'real'
    d.config.error_func = @error_abs;
    d.config.nonprojective_error = true;
    %d.config.UL_mixed = true;  % requires a task with mixed initial/final states

  case 'tr'
    d.config.error_func = @error_tr;

  case 'full'
    d.config.error_func = @error_full;

  otherwise
    disp('Keeping the old error function and gradient.')
    ttt = '';
end

% required after changing error/gradient func
d.cache_init()

% full gradient, including tau
mask = d.full_mask(true);


switch test
%% walltime benchmark
  case 'time'

% save the initial controls
%x0 = d.seq.get_raw(mask);
%d.update_controls(x0 + delta, mask);

tic
for k=1:10
    [err, grad] = d.compute_error(mask);
end
t = toc


%% test the gradient accuracy
  case 'acc'

% for flushing out gradient_setup bugs:
% perturb controls, do/do not recompute entire cache
%x1 = x0 +randn(size(x0));
%d.update_controls(x1, mask);
%d.cache_fill();  % recompute everything

% save the initial controls
x0 = d.seq.get_raw(mask);

% error function and its gradient at x0
[err, grad] = d.compute_error(mask);

if nargin < 3
    % random unit direction in parameter space
    direction = randn(size(x0));

elseif isempty(direction)
    % follow the gradient
    direction = grad;
end
direction = direction / norm(direction);
    
s = logspace(0, -10, 40);
diff = [];
predicted = [];
accurate = [];
for k=1:length(s)
    delta = s(k) * direction;

    % linear estimate of f(x) at x0+delta
    predicted(k) = err + grad.' * delta;
    
    % f(x0+delta)
    d.update_controls(x0 + delta, mask);
    accurate(k) = d.compute_error();
end

% restore initial controls
d.update_controls(x0, mask);


%% plot the results

diff = abs(predicted -accurate);

figure();
subplot(1,2,1)
% gradient error should be \propto s^2 for an exact or
% finite_diff gradient with a sufficiently small epsilon, and
% \propto s for a 1st order gradient approximation.
loglog(s, diff, 'b-o', s, s.^2, 'r', s, s, 'g');
xlabel('|\Delta x|')
ylabel('|error|')
grid on
legend('gradient error', 'quadratic', 'linear')
title(ttt)

subplot(1,2,2)
semilogx(s, predicted, 'b-o', s, accurate, 'r-o');
xlabel('|\Delta x|')
ylabel('f(x)')
grid on
legend('predicted', 'accurate');
title(ttt)
end
