% Example 5: S+E gate optimization

randseed(23483);

% 3-qubit Heisenberg chain, XY controls, QFT gate on the first 2 qubits
q = 3;
dim = 2 * ones(1, q);

desc = sprintf('Isotropic Heisenberg chain, %d qubits, XY control.', q);
fprintf('%s\n\n', desc);

J = 2 * [1 1 1]; % Heisenberg interaction
C = diag(ones(1, q-1), 1); % topology: linear chain

%C(2,3) = 0
H_drift = heisenberg(dim, @(s,a,b) J(s)*C(a,b));
[H_ctrl, c_labels] = control_ops(dim, 'xy');

final = qft(2);
initial = eye(prod(dim));

dyn = dynamo('closed gate_partial', initial, final, H_drift, H_ctrl);
dyn.system.set_labels(desc, dim, c_labels);
dyn.seq_init(20, 5 * [1, 0]);
dyn.set_controls(0.1);

dyn.ui_open();
dyn.search();
%dyn.analyze();
