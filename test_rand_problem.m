function dyn = test_rand_problem(task, dim, n_controls)
% Generates a Dynamo instance corresponding to a random optimization problem.
%
%  dyn = test_rand_problem(task, dim, n_controls)

% Ville Bergholm 2015


% S dimension
d = dim(1);

% SE total dimension
dSE = prod(dim);


switch task
  case {'closed ket', 'closed ket phase'}
    % pure state transfer
    ini = rand_ket(d);
    fin = rand_ket(d);

  case 'closed state'
    % mixed state transfer
    ini = rand_positive(d);
    fin = rand_positive(d);

  case {'closed gate', 'closed gate phase'}
    % unitary gate
    ini = eye(d);
    fin = rand_U(d);

  case 'closed state_partial'
    % mixed state transfer on SE, only S matters
    ini = rand_positive(dSE);
    fin = rand_positive(d);

  case 'closed gate_partial'
    % unitary gate on SE, only S matters
    ini = eye(dSE);
    fin = rand_U(d);

  case {'open state', 'open state overlap'}
    % mixed state transfer
    ini = rand_positive(d);
    fin = rand_ket(d);

  case 'open gate'
    % unitary gate
    ini = eye(d);
    fin = rand_U(d);

  case 'open state_partial'
    % mixed state transfer on SE, only S matters
    ini = rand_positive(dSE);
    fin = rand_positive(d);

  case 'open gate_partial'
    % TODO
    error('not finished')

  case 'abstract vector'
    ini = rand_complex([d,1]);
    fin = rand_complex([d,1]);

  case 'abstract matrix'
    ini = rand_complex(d);
    fin = rand_complex(d);

  otherwise
    error('Unknown task.')
end


H_ctrl = cell(1, n_controls);

switch strtok(task)
  case 'abstract'
    % abstract linear evolution
    H_drift = rand_complex(dSE);
    for k=1:n_controls
        H_ctrl{k} = rand_complex(dSE);
        c_labels{k} = sprintf('C%d', k);
    end

  case 'closed'
    % Hamiltonian evolution
    H_drift = rand_hermitian(dSE);
    for k=1:n_controls
        H_ctrl{k} = rand_hermitian(dSE);
        c_labels{k} = sprintf('H%d', k);
    end

  case 'open'
    % Markovian quantum evolution
    H_drift = rand_hermitian(dSE);
    % Lindblad generators
    for k=1:4
        A{k} = rand_complex(dSE) * 0.1;
    end
    H_drift = superop_lindblad(A, H_drift);
    for k=1:n_controls
        H_ctrl{k} = rand_hermitian(dSE);
        c_labels{k} = sprintf('H%d', k);
    end

  otherwise
    error('Unknown case.')
end

% control limits
%control_type = 'mm';
%temp = [-1,2] * 10;
%control_par = {temp, temp};




%% Set up Dynamo

desc = task;
T = 3;
n_bins = 100;

dyn = dynamo(task, ini, fin, H_drift, H_ctrl);
dyn.system.set_labels(desc, dim, c_labels);
dyn.seq_init(n_bins, T * [1, 0]); %, control_type, control_par);

% random, constant initial controls
dyn.easy_control(0.1 * randn(1, n_controls));

%dyn.ui_open();
%dyn.search(dyn.full_mask(false));
end


function ket = rand_ket(dim)
% Returns a random Haar-distributed ket vector.
    ket = zeros(dim, 1);
    ket(1) = 1;
    ket = rand_U(dim) * ket;
end

function ret = rand_complex(size)
% Returns a random complex matrix.
    ret = randn(size) +1i*randn(size);
end
