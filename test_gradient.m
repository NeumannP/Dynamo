function [d, direction] = test_gradient(d, direction)
% Checks if the computed gradient of the error function is accurate.
%
%  d is a Dynamo instance containing the optimization problem used.
%  If no d is given, uses one of the test suite problems.
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
    
% Ville Bergholm 2011-2015


%randseed(seed);


%% set up a system, random controls

if nargin < 1
    d = test_suite(21);
    d.easy_control({d.seq.fields, d.seq.tau+0.01*randn(size(d.seq.tau))});
end
   

%% choose an error function and a compatible gradient

ff = 'full'
gg = 'exact'

% for the finite_diff methods only
d.config.epsilon = 1e-3;

ttt = ['error\_', ff, ', gradient\_', gg];
switch ff
  case 'g'
    d.config.error_func = @error_abs;
    %d.config.error_func = @error_real;
    switch gg
      case 'exact'
        d.config.gradient_func = @gradient_g_exact;
      case '1st'
        d.config.gradient_func = @gradient_g_1st_order;
      case 'diff'
        d.config.gradient_func = @gradient_g_finite_diff;
      otherwise
        error('zzzz')
    end

  case 'tr'
    d.config.error_func = @error_tr;
    switch gg
      case 'exact'
        d.config.gradient_func = @gradient_tr_exact;
      case 'diff'
        d.config.gradient_func = @gradient_tr_finite_diff;
      otherwise
        error('zzzz')
    end

  case 'full'
    d.config.error_func = @error_full;
    switch gg
      case 'exact'
        d.config.gradient_func = @gradient_full_exact;
      case '1st'
        d.config.gradient_func = @gradient_full_1st_order;
      case 'diff'
        d.config.gradient_func = @gradient_full_finite_diff;
      otherwise
        error('zzzz')
    end

  otherwise
    disp('Keeping the old error function and gradient.')
    ttt = '';
end



%% test the accuracy

mask = d.full_mask(true);

% save the initial controls
x0 = d.seq.get(mask);

% error function and its gradient at x0
[err, grad] = d.compute_error(mask);

if nargin < 2
    % random unit direction in parameter space
    direction = randn(size(x0));

elseif isempty(direction)
    % follow the gradient
    direction = grad;
end
direction = direction / norm(direction);
    
s = logspace(0, -6, 30);
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
legend('gradient error', 'quadratic', 'linear')
title(ttt)

subplot(1,2,2)
semilogx(s, predicted, 'b-o', s, accurate, 'r-o');
xlabel('|\Delta x|')
ylabel('f(x)')
legend('predicted', 'accurate');
title(ttt)
end
