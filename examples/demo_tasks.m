function dyn = demo_tasks(task)
% Example: Demonstrate all the optimization tasks.
%
%  dyn = demo_tasks(task)

% Ville Bergholm 2014-2016


%% Pauli matrices

X = [0, 1; 1, 0];
Y = [0, -1i; 1i, 0];
Z = [1, 0; 0, -1];


%% Define the physics of the problem
% q-qubit Heisenberg chain

q = 2;                 % number of qubits
D_sub = 3;             % dimension of the interesting subspace (if applicable)
qs = 1;                % number of system qubits (if applicable)
dim = 2 * ones(1, q);  % dimension vector
D = prod(dim);         % total dimension

desc = sprintf('Isotropic Heisenberg chain, %d qubits, XY control at one end, task: %s', q, task);

% Drift Hamiltonian
J = 2 * pi * [1 1 1]; % Heisenberg interaction
C = diag(ones(1, q-1), 1); % topology: linear chain
H_drift = heisenberg(dim, @(s,a,b) J(s)*C(a,b));

%H_drift = H_drift -op_sum(dim, @(k) (k+2)*Z);


% XY controls at one end of the chain
[H_ctrl, c_labels] = control_ops(dim, '1xy');

%H_ctrl = {op_sum(dim, X), op_sum(dim, Y)};
%c_labels = {'X', 'Y'};

% limit driving Rabi frequency to the interval [-1, 1]
control_type = '..';
temp = [-1,2] * 10;
control_par = {temp, temp};


% adjust the physics depending on the task
switch strtok(task)
  case 'abstract'
    % our abstract system example uses a non-hermitian generator
    dim = 3;
    desc = sprintf('Non-hermitian three-level system, task: %s', task);

    % abstract control operators (not Hamiltonians!)
    C1 = sparse(3,3);
    C1(1,2) = -1;
    C1(2,1) = 1;
    C2 = sparse(3,3);
    C2(2,3) = -1;
    C2(3,2) = 1;
    H_ctrl = {C1, C2};
    c_labels = {'Omega_p', 'Omega_s'};

    % drift operator
    H_drift = sparse(3,3);
    H_drift(2,2) = -1;

  case 'closed'

  case 'open'
    % in an open system, add some dephasing on each qubit
    temp = control_ops(dim, 'd');
    L_deph = 0;
    for k=1:q
        L_deph = L_deph +0.1 * temp{k};
    end
    % change the drift Hamiltonian into a Liouvillian
    H_drift = superop_lindblad({}, H_drift) +L_deph;
  
  otherwise
    error('Unknown task.')
end


%% Choose the task
% The initial and final states depend on the optimization task.

initial = zeros(D,1);  initial(1) = 1;
final = zeros(D,1);  final(end) = 1;

switch task
  case {'closed ket', 'closed ket phase'}
    % pure state transfer
    T = 0.72;

  case 'closed state'
    % mixed state transfer
    a = 0.6;
    initial = a * initial * initial' +(1-a) * eye(D) / D;
    T = 0.72;
    
  case {'closed gate', 'closed gate phase'}
    % QFT gate
    initial = eye(D);
    final = qft(q);
    T = 1.6;

  case 'closed gate subspace'
    % random gate on a subspace
    initial = eye(D);
    % pad the non-interesting parts of final with zeroes
    P = [eye(D_sub), zeros(D_sub, D-D_sub)];
    final = P' * rand_U(D_sub) * P
    T = 3;

  case 'closed state_partial'
    % first qubit to |1>
    final = [0; 1];
    T = 0.59;

  case 'closed gate_partial'
    % partial QFT gate
    initial = eye(D);
    final = qft(qs);
    T = 1.3;

  case {'open state', 'open state overlap'}
    T = 1;
    
  case 'open gate'
    initial = eye(D);
    final = qft(q);
    T = 1.8;

  case 'open state_partial'
    % first qubit to |1>
    final = [0; 1];
    T = 0.59;

  case 'open gate_partial'
    % partial QFT gate
    initial = eye(D);
    final = qft(qs);
    T = 1.3;

  case 'abstract vector'
    initial = [1 0 0].';
    final = [0 0 1].';
    T = 0.7;
    
  case 'abstract matrix'
    initial = eye(3);
    final = [1, 0, 0; 0, 0, 1; 0, -1, 0] / sqrt(3); % TODO kinda arbitrary
    T = 1.5;

  otherwise
    error('Unknown task.')
end


%% Set up Dynamo

fprintf('%s\n\n', desc);
bins = 100;

dyn = dynamo(task, initial, final, H_drift, H_ctrl);
dyn.system.set_labels(desc, dim, c_labels);
dyn.seq_init(bins, T * [1, 0], control_type, control_par);

% random, constant initial controls
dyn.shake(0, 0.4, false);

% do not optimize taus
mask = dyn.full_mask(false);


%% Now do the actual search

dyn.ui_open();
dyn.search(mask);
%dyn.analyze();
%figure; dyn.plot_pop();
end
