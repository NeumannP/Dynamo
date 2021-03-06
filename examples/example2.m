% Example 2 for DPG 2012 Stuttgart

randseed(2783);

% 3-qubit Heisenberg chain, XY controls at one end, QFT gate
q = 3;
dim = 2 * ones(1, q);

desc = sprintf('Isotropic Heisenberg chain, %d qubits, XY control at one end.', q);
fprintf('%s\n\n', desc);

J = 2 * [1 1 1]; % Heisenberg interaction
C = diag(ones(1, q-1), 1); % topology: linear chain
H_drift = heisenberg(dim, @(s,a,b) J(s)*C(a,b));
[H_ctrl, c_labels] = control_ops(dim, '1xy');

final = qft(q);
initial = eye(size(final));

dyn = dynamo('closed gate', initial, final, H_drift, H_ctrl);
dyn.system.set_labels(desc, dim, c_labels);
dyn.seq_init(100, 16 * [1, 0]);
dyn.set_controls(0.1);

dyn.ui_open();
dyn.search();
%dyn.analyze();
