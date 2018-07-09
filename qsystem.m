classdef qsystem < matlab.mixin.Copyable
% Copyable handle class defining a quantum system.

% Ville Bergholm 2011-2016
    
  properties
    description = ''    % description string
    dim                 % dimension (or dim vector) of the Hilbert space of the full system S+E
    dimSE               % total dimensions of S (the part we're interested in) and the environment E, as a two-vector
    type = 'hilbert';   % Do the system objects X reside in Liouville or Hilbert space?

    weight = 1          % vector, length == n_ensemble, weights for the ensemble samples
    A                   % cell vector of drift generators, size == [n_ensemble, 1]
    B                   % cell array of control generators, size == [n_ensemble, n_controls]
    B_is_Hamiltonian    % vector, size == [1, n_controls], true when corresponding B items represent Hamiltonians

    X_initial           % initial state
    X_final             % final state
    norm2               % squared norm of final state
    TU = []             % time unit for generators, in seconds. \hat{A} = A * TU etc.
    state_labels   = {} % names for the Hilbert space computational basis states
    control_labels = {} % names for the controls
  end


  methods (Static)
      function ret = is_hermitian(H)
      % Returns true iff H is hermitian
          ret = (norm(H -H') < 1e-10);  % FIXME magic number
      end

      function ret = get_A(A, k)
      % Returns the drift generator A corresponding to ensemble index k.
      % This and the corresponding get_B function only exist to
      % enable more flexible input of A and B.

          if isa(A, 'function_handle')
              ret = A(k);
          elseif iscell(A)
              % cell array of matrices
              ret = A{k};
          else
              % single matrix
              ret = A;
          end
          % eig and norm cannot handle sparse matrices...
          ret = full(ret);
      end

      function ret = get_B(B, k, c, n_ensemble)
      % Returns the control generator B corresponding to ensemble index k, control field c.

          if isa(B, 'function_handle')
              ret = B(k, c);
          else
              % cell array of matrices
              ret = B{k, c};
          end
          ret = full(ret);
      end
  end


  methods (Access = private)
    function [ret, is_H] = setup_abstract(self, C, message)
    % Process abstract generators.
        ret = C;
        is_H = false;
    end

    function [ret, is_H] = setup_hamiltonian(self, H, message)
    % Process Hamiltonian generators. Raises an error if H is not a valid Hamiltonian.
        if ~qsystem.is_hermitian(H)
            error(strcat(message, ' is not hermitian.'))
        end
        ret = -1i * H;
        is_H = true;
    end

    function [ret, is_H] = setup_liouvillian(self, C, message)
    % Convert Hamiltonians into Liouvillians if necessary.

        if length(C) == self.dim
            % assume it's a Hamiltonian
            [ret, is_H] = self.setup_hamiltonian(C, message);
            ret = full(superop_comm(ret));  % into a commutator superoperator
        else
            % assume it's a Liouvillian
            ret = C;
            is_H = false;
        end
    end

    function setup_gens(self, A, B, func)
    % Check and set up drift and control generators.

        if iscell(A) && length(A) ~= self.n_ensemble()
            error('Number of elements in cell array A does not match n_ensemble.')
        end
        if iscell(B)
            s = size(B);
            if s(1) == 1
                B = @(k,c) B{1,c};  % same for every ensemble member
            elseif s(1) ~= self.n_ensemble()
                error('Number of rows in cell array B does not match n_ensemble.')
            end
            if s(2) ~= self.n_controls()
                error('Specified number of controls does not match the number of columns in cell array B.')
            end
        elseif ~isa(B, 'function_handle')
            error('Control generators must be given as a cell array or as a function handle.');
        end

        % Store the generators internally as cell arrays for efficiency.
        temp = true(size(self.B));
        for k = 1:self.n_ensemble()
            msg = sprintf('Drift Hamiltonian (ensemble %d)', k);
            self.A{k} = func(self, qsystem.get_A(A, k), msg);
            for c = 1:self.n_controls()
                msg = sprintf('Control Hamiltonian %d (ensemble %d)', c, k);
                [self.B{k, c}, temp(k, c)] = func(self, qsystem.get_B(B, k, c), msg);
            end
        end
        % TODO it does not make much sense to allow the type of a single
        % control to vary between ensemble members... we should maybe check it here
        self.B_is_Hamiltonian = all(temp, 1);
    end
  end
  
  
  methods
    function self = qsystem(weight, input_dim, n_controls)
    % constructor
        
        % ensemble init
        if ~isvector(weight)
            error('Ensemble weights must be given in a vector.')
        end
        self.weight = weight;
        n_ensemble = self.n_ensemble();

        % Set up the Hilbert space dimensions.
        % E may not exist, in which case it has dimension 1.
        self.dim      = input_dim(1); % S+E, dim of initial
        self.dimSE(1) = input_dim(2); % S, dim of final
        temp = self.dim / self.dimSE(1); % E
        if floor(temp) ~= temp  % must be an integer
            error('Initial state must be an object on S+E, final state an object on S.');
        end
        self.dimSE(2) = temp;

        % initialize the generator containers
        self.A = cell(n_ensemble, 1);
        self.B = cell(n_ensemble, n_controls);
        self.B_is_Hamiltonian = true(1, n_controls);
    end


    function abstract_representation(self, i, f, A, B)
    % X_ are abstract Hilbert space objects (vectors or matrices).

        self.type = 'abstract';
        self.X_initial = i;
        self.X_final = f;
        self.norm2 = norm2(self.X_final);

        % set the generators
        self.setup_gens(A, B, @setup_abstract)
    end


    function hilbert_representation(self, i, f, A, B, gate_partial)
    % X_ are Hilbert space objects (kets or operators).
    % Used for closed system tasks: state, ket, gate, gate_partial, (TODO state_partial).
    % For _partial tasks, i \in SE, f \in S.
    % (NOTE: the generators are not pure Hamiltonians, there's an extra -1i!)
        
        self.type = 'hilbert';
        self.X_initial = i;
        if gate_partial
            % only with gate_partial
            self.X_final = kron(f, eye(self.dimSE(2)));
        else
            self.X_final = f;
        end
        
        % Calculate the squared norm |X_final|^2 to scale the fidelities with.
        % We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
        self.norm2 = norm2(self.X_final);
    
        % set up the generators
        self.setup_gens(A, B, @setup_hamiltonian);
    end


    function vec_representation(self, i, f, A, B, use_states)
    % X_ are Liouville space vectors/operators corresponding to vec-torized state operators / unitary gates.
    % Used for open system tasks: state, state_partial, gate, (TODO gate_partial).

        self.type = 'liouville';
        if use_states
            % state vectors are converted to state operators
            self.X_initial = vec(to_op(i));
            self.X_final   = vec(to_op(f));
        else
            % i and f are gates
            self.X_initial = lrmul(i, i'); % == kron(conj(i), i);
            self.X_final   = lrmul(f, f'); % == kron(conj(f), f);
        end
        self.norm2 = norm2(self.X_final);

        % Set up Liouville space generators.
        self.setup_gens(A, B, @setup_liouvillian);
    end        


    function set_TU(self, TU)
    % Sets the time unit for the system.
        self.TU = TU;
    end


    function set_labels(self, desc, st_labels, c_labels)
    % Describe the system, label the states and controls. The labels are cell vectors of strings.

        self.description = desc;
        D = self.dimSE(1); % label just the S states
        n_controls = self.n_controls();
        
        if nargin < 3 || isempty(st_labels)
            % use default state labels
            st_labels = char('0' + (1:D).');
        elseif ~iscell(st_labels)
            % it's a dim vector, use standard computational basis labeling
            dim = st_labels;
            if prod(dim) ~= prod(self.dimSE)
                error('Dimension vector given does not match the Hilbert space dimension.')
            end
            self.dim = dim; % store the dim vector (replacing the default scalar total dimension)

            % find where S ends and E starts
            n = find(cumprod(dim) == D, 1);
            ket = zeros(1,n);
            st_labels = cell(1, D);
            % build the labels
            for k=1:D
                st_labels{k} = ['$|', char(ket + '0'), '\rangle$'];
                for b = n:-1:1 % start from least significant digit
                    ket(b) = ket(b)+1;
                    if ket(b) < dim(b)
                        break;
                    end
                    ket(b) = 0;
                end
            end
        end

        if nargin < 4 || isempty(c_labels)
            c_labels = char('0' + (1:n_controls).');
        end
        
        if length(st_labels) ~= D
            error('Number of state labels given does not match the Hilbert space dimension of S.')
        end
        self.state_labels = st_labels;
        
        if length(c_labels) ~= n_controls
            error('Number of control labels given does not match the number of controls.')
        end
        self.control_labels = c_labels;
    end
    
    
    function ret = n_ensemble(self)
    % Returns the number of systems in the ensemble sample.
        ret = length(self.weight);
    end

    function ret = n_controls(self)
    % Returns the number of control fields.
        ret = size(self.B, 2);
    end
  end
end
