function dyn = demo_abstract(T)
% Example: State transfer for an abstract three-level system.

% Robert Zeier 2013
% Ville Bergholm 2013-2014


%% define the problem

dim = 3;
desc = sprintf('State transfer for a three-level system.');
fprintf('%s\n\n', desc);

% abstract control operators (not Hamiltonians!)
C1 = sparse(3,3);
C1(1,2) = -1;
C1(2,1) = 1;
C2 = sparse(3,3);
C2(2,3) = -1;
C2(3,2) = 1;

c_labels = {'Omega_p', 'Omega_s'};

% drift operator
H_drift = sparse(3,3);
H_drift(2,2) = -1;

final = [0 0 1].';
initial = [1 0 0].';

dyn = dynamo('abstract vector', initial, final, H_drift, {C1, C2});
dyn.system.set_labels(desc, dim, c_labels);

dyn.config.error_func = @error_full;
dyn.config.gradient_func = @gradient_full_finite_diff;
dyn.config.epsilon = 1e-4;


%% initial controls

if nargin < 1
    T = 0.5;
end

% control sequence consists of n_b = 100 bins (first parameter),
% the time duration (tau) of each bin fixed to 1/n_b time units (second parameter).
dyn.seq_init(100, T * [1, 0]);
% The second parameter defines the tau limits for each bin: [tau_minimum, tau_delta]
% so that tau is always between tau_minimum and tau_minimum + tau_delta.
% You can give an [n_b, 2] -sized array (separate limits for each bin),
% or just a [1,2]-sized vector (like here) which is then divided by n_b and applied to each bin.
% The tau values are by default initialized to the middle of their respective intervals.
% In this case tau_delta is zero so the taus are fixed.


% Initialize the controls (to constant, random values).
% These aren't limits for the controls, just the initial values.
%dyn.easy_control([-10, 10]);
dyn.easy_control(randn(1,2));

% Control limits (optional) can be given in the seq_init call, like this:
%dyn.seq_init(n_b, T * [0.5, 1], c_type, c_par);
% See control_seq.m:set()


%% optimize

dyn.ui_open();
dyn.search();
%dyn.analyze();
%figure; dyn.plot_X();
%figure; dyn.plot_seq();
