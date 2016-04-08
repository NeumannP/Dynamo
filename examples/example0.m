% Example 0 for DPG 2012 Stuttgart

randseed(2783);

% 2-qubit Heisenberg chain, XY controls at one end, state transfer |00> -> |11>

q = 2;
dim = 2 * ones(1, q);

desc = sprintf('Isotropic Heisenberg chain, %d qubits, XY control at one end.', q);
fprintf('%s\n\n', desc);

J = 2 * [1 1 1];
C = diag(ones(1, q-1), 1); % topology: linear chain
H_drift = heisenberg(dim, @(s,a,b) J(s)*C(a,b));
[H_ctrl, c_labels] = control_ops(dim, '1xy');

initial = [1 0 0 0].';
final = [0 0 0 1].';

dyn = dynamo('closed ket', initial, final, H_drift, H_ctrl);
dyn.system.set_labels(desc, dim, c_labels);
dyn.seq_init(100, 6 * [1, 0]);
dyn.easy_control([-0.1, 0.05]);

dyn.ui_open();
dyn.search();
%dyn.analyze();
%figure; dyn.plot_X();
%figure; dyn.plot_seq();
