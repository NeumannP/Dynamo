function dyn = demo_resonance(delta)
% Example: DYNAMO discovers the importance of resonant driving on its own.
%
%  dyn = demo_resonance(delta)
%
%  Optimizes a single-qubit control sequence for the state transfer |0> \to |1>.
%  delta is the detuning (difference between the carrier and
%  resonance frequencies). Uses RWA, ignores the fast-rotating term.
%
%  Resonant driving (delta = 0) yields just a constant-strength pi pulse.
%  Slightly off-resonant driving (delta = 0.1) discovers a CORPSE-type sequence.
%  Far off-resonant driving (delta = 5) results in a sinusoidal pulse that
%  tries to approximate resonant driving.

% Ville Bergholm 2011-2014

if nargin < 1
    delta = 5;
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
dyn.system.set_labels('Single-qubit resonant driving demo.', dim, c_labels);


%% Initial controls

% random initial controls
T = 6;
dyn.seq_init(201, T * [1, 0], control_type, control_par);
dyn.easy_control(0, 0, 1, false);


%% Now do the actual search

dyn.ui_open();
dyn.search();
dyn.analyze();
end
