function dyn = demo_tau(delta)
% Example: Optimizing the bin durations.
%
%  dyn = demo_tau(delta)
%
%  Optimizes a single-qubit control sequence for the state transfer |0> \to |1>.
%  The optimizer can compensate for the very limited number of bins by
%  changing their durations.

% Ville Bergholm 2014-2015

if nargin < 1
    delta = 0.6;
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
[H_ctrl, c_labels] = control_ops(dim, 'x');

% limit driving Rabi frequency to the interval [-1, 1]
control_type = 'm';
control_par = {[-1,2]};

% pure state transfer
initial = [1 0]';
final = [0 1]';

dyn = dynamo('closed ket', initial, final, H_drift, H_ctrl);
dyn.system.set_labels('Bin duration optimization demo.', dim, c_labels);


%% Initial controls

% going for a CORPSE-style sequence
T = 7/3 * pi;
% allow taus to vary between 0.5 and 1.5 times the original duration
dyn.seq_init(3, T * [0.5, 1], control_type, control_par);
% random initial controls
dyn.set_controls(rand());

% optimize taus as well
mask = dyn.full_mask(true);
%mask(:,1) = 0


%% Now do the actual search

dyn.ui_open();
dyn.search(mask, 'Display', 'final', 'plot_interval', 1);
%dyn.analyze();
end
