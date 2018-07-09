function dyn = demo_ensemble(n_ensemble, delta)
% Example: Ensemble optimization.
%
%  dyn = demo_ensemble(n_ensemble, delta)
%
%  Optimizes a robust single-qubit control sequence for the state transfer |0> \to |1>
%  over an ensemble of systems with slightly different detunings.
%  Uses RWA, ignores the fast-rotating term.

% Ville Bergholm 2014-2016

if nargin < 2
    delta = 1;
    if nargin < 1
        n_ensemble = 5;
    end
end

%% Pauli matrices etc.

Z = diag([1, -1]);


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

% Drift Hamiltonian
H_drift = -delta * Z / 2;

% Control Hamiltonians / Liouvillians
[H_ctrl, c_labels] = control_ops(dim, 'xy');
n_controls = length(H_ctrl);

% limit driving Rabi frequency to the interval [-1, 1]
control_type = 'mm';
pp = [-1, 2];
control_par = {pp, pp};

% pure state transfer
initial = [1 0]';
final = [0 1]';

% define the ensemble
if 1
    str = 'detuning';
    x = linspace(-1, 1, n_ensemble);
    sigma = 0.5;
    A_func = @(x)   x * H_drift;  % vary drift
    B_func = @(x,c) H_ctrl{c};
else
    str = 'relative pulse strength error';
    x = linspace(-0.1, 0.1, n_ensemble)
    sigma = 0.08;
    A_func = @(x)   0.1 * H_drift;
    B_func = @(x,c) (1+x) * H_ctrl{c}; % vary control scaling
end
% gaussian distribution of weights
weight = exp(-x.^2 / (2*sigma^2));
weight = weight / sum(weight)

dyn = dynamo('closed ket', initial, final, @(k) A_func(x(k)), @(k,c) B_func(x(k),c), weight, n_controls);
dyn.system.set_labels('Single-qubit ensemble optimization demo.', dim, c_labels);


%% Initial controls

% random initial controls
T = pi * 7/3;  % enough for the short CORPSE sequence
dyn.seq_init(201, T * [1, 0], control_type, control_par);
dyn.set_controls(2*rand(1, 2)-1);


%% Now do the actual search

dyn.ui_open();
dyn.search();
%dyn.analyze();


%% Plot the sequence error over different ensemble parameter values

xxx = 2*linspace(x(1), x(end), 100);
err = [];
for k=1:length(xxx)
    % keep the sequence, get rid of the ensemble, change the generators
    dyn.system = qsystem(1, [size(initial, 1), size(final, 1)], n_controls);
    dyn.system.hilbert_representation(initial, final, A_func(xxx(k)), @(dummy,c) B_func(xxx(k),c), false);
    dyn.cache_init()
    err(k) = dyn.compute_error();
end
err = sqrt(2 * dyn.system.norm2 * err); % into Frobenius norm error
figure
title('Sequence error as a function of the ensemble parameter')
xlabel(str);
ylabel('Frobenius norm error');
axis([xxx(1), xxx(end), 0, max(err)]);
grid on;
hold on;
plot(xxx, err);
plot(x, min(err), 'ko')
legend('sequence error', 'ensemble')
end
