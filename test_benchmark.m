function [d] = test_benchmark(d)
% Measures the walltime required to compute error functions and gradients.
%
%  d is a Dynamo instance containing the optimization problem used.
%  If no d is given, uses a random problem.

% Ville Bergholm 2015


%% set up a system, random controls

if nargin < 1
    d = test_rand_problem('closed gate', 4, 2);
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



%% benchmark

mask = d.full_mask(true);

% save the initial controls
%x0 = d.seq.get(mask);
%d.update_controls(x0 + delta, mask);

tic
for k=1:10
    [err, grad] = d.compute_error(mask);
end
toc
end
