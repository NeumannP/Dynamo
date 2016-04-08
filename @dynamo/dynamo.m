classdef dynamo < matlab.mixin.Copyable
% Copyable handle class for DYNAMO optimizer objects.
%
% Contains the optimization task, system description, control
% sequence, various options and statistics. This class glues together
% the functionality of the cache, qsystem and control classes.
%
% Governing equation: \dot(X)(t) = (A +\sum_k u_k(t) B_k) X(t) = G(t) X(t)
    
% Shai Machnes   2010-2011
% Ville Bergholm 2011-2015


  properties
    config   % configuration information, other metadata
    system   % description of the physical system
    seq      % control sequence
    opt      % optimization options
    stats    % optimization statistics
  end

  properties (Transient)
    cache    % Do not save the cache on disk since it may be huge and can always be recomputed.
  end

  methods (Static)
    function ret = version()
    % Returns the current DYNAMO version.
        ret = '1.4.0alpha4';
    end


    function obj = loadobj(obj)
    % Re-initializes the cache (which is not saved) during loading.
        if isstruct(obj)
            error('Backwards compatibility of saved objects not yet implemented.')
        end
        % HACK, backwards compatibility. Assume that dim(E) = 1.
        if isempty(obj.system.dimSE)
            obj.system.dimSE = [prod(obj.system.dim), 1];
        end
        obj.cache_init();
    end
  end

  methods (Access = protected)
    function cp = copyElement(self)
    % Override the default copyElement method to provide deep copies.
    % This is necessary since some of the data members are Copyable handle classes themselves.
        
        % Make a shallow copy of everything
        cp = copyElement@matlab.mixin.Copyable(self);
        % Make a deep copy of all handle-type properties
        cp.system = copy(self.system);
        cp.seq    = copy(self.seq);
        cp.cache  = copy(self.cache);
    end
  end
  
  methods
    function self = dynamo(task, initial, final, A, B, weight)
    % Constructor

        if nargin < 6
            % ensemble weights
            weight = 1;
        end
        
        task = lower(task);

        %% Some basic data provenance

        config.version = dynamo.version();
        % Local time. UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
        config.date = datestr(now(), 31);
        % HACK: UTC time using Java
        temp = datenum('1970', 'yyyy') +java.lang.System.currentTimeMillis / 86400e3;
        config.date_UTC = datestr(temp, 31);
        
        config.task = task;
        config.UL_hack = false;
        
        [system_str, rem] = strtok(task);
        [task_str, rem] = strtok(rem);
        [extra_str, rem] = strtok(rem);
          

        %% Description of the physical system

        input_dim  = [size(initial, 1), size(final, 1)];
        input_rank = [size(initial, 2), size(final, 2)]; % check the validity of the inputs
        sys = qsystem(weight, input_dim, size(B, 2)); % FIXME assumes B is a cell array
        
        out = 'Target operation:';
        switch system_str
          case 'abstract'
            %% No transformations done on the A and B operators. 
            % the generator may be anything, hence error_full
            out = strcat(out, ' abstract');
            
            if strcmp(task_str, 'vector')
                out = strcat(out, ' vector transfer\n');
                if any(input_rank ~= 1)
                    error('Initial and final states should be vectors.')
                end
            else
                out = strcat(out, ' matrix operation\n');
                if any(input_rank == 1)
                    error('Initial and final states should be matrices.')
                end
            end
            sys.abstract_representation(initial, final, A, B);
            config.error_func = @error_full;
            config.gradient_func = @gradient_full_finite_diff;
            config.epsilon = 1e-4;
            
          
          case {'closed'}
            %% Closed system
            % the generator is always Hermitian and thus normal => use exact gradient
            switch task_str
              case 'state'
                % TEST more efficient Hilbert space implementation
                out = strcat(out, ' mixed state transfer');
                if all(input_rank == 1)
                    error('Either the initial or the final state should be a state operator.')
                end
                sys.hilbert_representation(to_op(initial), to_op(final), A, B, false);
                config.f_max = (sys.norm2 +norm2(sys.X_initial)) / 2;
                config.error_func = @error_real;
                config.gradient_func = @gradient_g_mixed_exact;
                config.UL_hack = true;
                
              
              case {'ket', 'gate'}
                if strcmp(task_str, 'ket')
                    out = strcat(out, ' pure state transfer');
                    if any(input_rank ~= 1)
                        error('Initial and final states should be normalized kets.')
                    end
                else
                    out = strcat(out, ' unitary gate');
                    if any(input_rank == 1)
                        error('Initial and final states should be unitary operators.')
                    end
                end
                sys.hilbert_representation(initial, final, A, B, false);
                config.f_max = sys.norm2;

                if strcmp(extra_str, 'phase')
                    out = strcat(out, ' (with global phase (NOTE: unphysical!))');
                    config.error_func = @error_real;
                else
                    out = strcat(out, ' (ignoring global phase)');
                    config.error_func = @error_abs;
                end
                config.gradient_func = @gradient_g_exact;

                
              % system S + environment E
              case 'gate_partial'
                out = strcat(out, ' partial unitary gate (on S)');
                if any(input_rank == 1)
                    error('Initial and final states should be unitary operators.')
                end
                sys.hilbert_representation(initial, final, A, B, true);
                config.f_max = sys.norm2;
                config.error_func = @error_tr;
                config.gradient_func = @gradient_tr_exact;
                
              case 'state_partial'
                error('Not implemented yet.')
                
              otherwise
                error('Unknown task.')
            end
            
            out = strcat(out, ' in a closed system.\n');

            
          case {'open'}
            %% Open system with a Markovian bath
            % The generator isn't usually normal, so we cannot use the exact gradient method

            switch task_str
              case 'state'
                out = strcat(out, ' state transfer');
                sys.vec_representation(initial, final, A, B, true);
                if strcmp(extra_str, 'overlap')
                    % overlap error function
                    % NOTE simpler error function and gradient, but final state needs to be pure
                    out = strcat(out, ' (overlap)');
                    config.error_func = @error_real;
                    config.gradient_func = @gradient_g_1st_order;
                    config.f_max = sys.norm2;
                else
                    % full distance error function
                    config.error_func = @error_full;
                    config.gradient_func = @gradient_full_1st_order;
                end
                
              case 'gate'
                out = strcat(out, ' quantum gate');
                if any(input_rank == 1)
                    error('Initial and final states should be unitary operators.')
                end
                sys.vec_representation(initial, final, A, B, false);
                config.error_func = @error_full;
                config.gradient_func = @gradient_full_1st_order;


              % system S + environment E
              case 'state_partial'
                out = strcat(out, ' partial state transfer (on S)');
                sys.vec_representation(initial, final, A, B, true);
                config.error_func = @error_full;
                config.gradient_func = @gradient_full_1st_order;
                
              case 'gate_partial'
                error('Not implemented yet.')

              otherwise
                % TODO arbitrary quantum maps
                error('Unknown task.')
            end
            
            out = strcat(out, ' in an open system under Markovian noise.\n');


          otherwise
            error('Unknown system specification.')
        end
        fprintf(out);
        if sys.liouville
            fprintf('Liouville space dimension: %d\n', sys.dim^2);
        else
            fprintf('Hilbert space dimension: %d\n', sys.dim);
        end
        n_ensemble = sys.n_ensemble();
        if n_ensemble > 1
            fprintf('Optimizing over an ensemble of %d systems.\n', n_ensemble);
        end
        fprintf('\n');
        
          
        % store the prepared fields
        self.config = config;
        self.system = sys;

        % init miscellaneous things
        self.config.ui_fig = [];
        self.stats = {};
    end


    function cache_init(self)
    % Set up cache after the number of time slots changes.
    % This is where all the bad code went.

        U_start = self.system.X_initial;

        % some error functions need a full reverse propagator.
        temp = self.config.error_func;
        if isequal(temp, @error_full)
            % NOTE X_initial, because in state_partial tasks X_final only describes the first subsystem
            L_end = eye(length(self.system.X_initial)); % L: full reverse propagator
        else
            L_end = self.system.X_final'; % L: X_final' propagated backwards
        end

        % exact gradient? we need the eigendecomposition data.
        temp = self.config.gradient_func;
        use_eig = isequal(temp, @gradient_g_exact)...
                  || isequal(temp, @gradient_g_mixed_exact)...
                  || isequal(temp, @gradient_tr_exact);

        % UL_hack: mixed states in a closed system
        self.cache = cache(self.seq.n_timeslots(), self.system.n_ensemble(), U_start, L_end, use_eig, self.config.UL_hack);
    end


    function seq_init(self, n_timeslots, tau_par, varargin)
    % Create the control sequence and a matching cache.
    % The varargin are the control_type and control_par cell vectors.

        self.seq = control_seq(n_timeslots, self.system.n_controls(), tau_par, varargin{:});
        if any(self.seq.control_type(~self.system.B_is_Hamiltonian) == '.')
            disp('Warning: Liouvillian control ops with possibly negative control values.')
        end
    
        self.cache_init();
    end


    function mask = full_mask(self, optimize_tau)
    % Returns a full control mask.
        
        n_timeslots = self.seq.n_timeslots();
        n_controls = self.seq.n_controls();

        % shape vectors
        f_shape = [n_timeslots, n_controls];
        t_shape = [n_timeslots, 1];

        %% Build the control mask

        fprintf('Tau values ');
        if optimize_tau
            fprintf('optimized.\n');
            mask = [true(f_shape), true(t_shape)];
        else
            fprintf('fixed.\n')
            mask = [true(f_shape), false(t_shape)];
        end
    end


    function [err_out, grad_out] = compute_error(self, control_mask)
    % Returns the error (and its gradient) at current control values.
    % This is where we sum over the system ensemble if necessary.

        % gradient requires a control mask
        if nargout == 2 && nargin < 2
            control_mask = self.full_mask(false);
        end

        if nargout == 2
            % since g can be computed using any U and L, it might be
            % cheaper to set up the gradient first...
            self.gradient_setup(control_mask);
        end

        % set up stuff for the error functions
        if isequal(self.config.error_func, @error_full)
            % _full:
            self.cache.g_needed_now = 2; % HACK ln37ae983e
        else
            % _tr, _abs, _real:
            self.cache.g_needed_now = true;
        end
        
        self.cache_refresh(); % this call does the heavy computing (expm etc.)

        err_out  = 0;
        grad_out = 0;
        % loop over the ensemble
        n_ensemble = self.system.n_ensemble();
        for k = 1:n_ensemble
            %% (real) normalized error
            err = self.config.error_func(self, k) / self.system.norm2;
            fprintf('Error (%d): %g\n', k, err);

            % weighted ensemble average
            err_out = err_out +self.system.weight(k) * err;
            if nargout < 2
                % just the error
                continue
            end

            %% gradient

            % tau are the last column of the controls
            tau_c = size(control_mask, 2);

            % iterate over the true elements of the mask
            grad = zeros(nnz(control_mask), 1);
            [Ts, Cs] = ind2sub(size(control_mask), find(control_mask));
            for z = 1:length(Ts)
                %fprintf('.');
                t = Ts(z);
                c = Cs(z);
                if c == tau_c
                    c = -1; % denotes a tau control
                end
                grad(z) = self.config.gradient_func(self, t, k, c);
            end
            % real, normalized, weighted gradient
            grad_out = grad_out +(self.system.weight(k) / self.system.norm2) * real(grad);
        end
        %fprintf('Error: %g\n', err_out);
    end

    
    function update_controls(self, raw, control_mask)
    % Updates selected controls.
    %
    %  raw: vector of raw, untransformed control values, control_mask: corresponding mask.
    %
    %  Updates control values for which control_mask is true.
    %  Makes the changed timeslots and stuff that depends on them stale.

        old = self.seq.get();
         
        if nargin < 3 || isempty(control_mask)
            control_mask = true(size(old)); % full mask
        end

        % make a trial copy of the new controls
        new = old;
        new(control_mask) = raw;

        % see which timeslots have changed
        changed_t_mask = any(new ~= old, 2);

        if any(changed_t_mask)
            % actually update the controls
            self.seq.set(new);
            self.cache.mark_as_stale(changed_t_mask);
        end
    end

    
    function cache_refresh(self)
    % Performs all the queued computations using the cache subsystem.
        self.cache.refresh(self.system, self.seq.tau, self.seq.fields);
    end
    

    function cache_fill(self)
    % Will invalidate everything, then re-calc everything in the cache.
    % Used mostly for debugging (since it essentially overrides all matrix-op optimization mechanisms).

        self.cache.invalidate();

        self.cache.H_needed_now(:) = true;
        self.cache.P_needed_now(:) = true;
        self.cache.U_needed_now(:) = true;
        self.cache.L_needed_now(:) = true;
        self.cache.g_needed_now    = true;
        
        self.cache_refresh();
    end
    
  
    function ret = X(self, j, k)
    % Returns X(t_j), the controlled system at time t_j.
    % If no j is given, returns the final state X(t_n).
    % k is an optional ensemble index.

        n_timeslots = self.seq.n_timeslots();

        if nargin < 3 || isempty(k)
            k = 1; % TODO for the lack of a better option. maybe average over the ensemble?
        end
        if nargin < 2 || isempty(j)
            j = n_timeslots; % final time
        end

        if k <= 0 || k > self.system.n_ensemble()
            error('Bad ensemble index.')
        end
        if j < 0 || j > n_timeslots
            error('Bad timeslot.')
        end
        % U{j+1} is the system at t = sum(tau(1:j)) = t_j
        self.cache.U_needed_now(j+1) = true;
        self.cache_refresh();
        ret = self.cache.U{j+1, k};
    end


    function plot_seq(self, varargin)
    % Plots the control sequence. Wrapper.

        ax = self.seq.plot(self.system.control_labels, varargin{:});
        title(ax, self.system.description);
        if ~isempty(self.system.TU)
            xlabel(ax, sprintf('time (%g s)', self.system.TU));
        end
    end


    function plot_stats(self, stat, ax)
    % Plots the optimization stats (like the error) as a function of wall time.

        fp = @plot;
        switch stat
          case 'error'
            fp = @(a,t,x,m) semilogy(a, t, abs(x), m);
          case 'control_integral'
          otherwise
            error('Unknown stat.')
        end

        offset = 0;
        for k=1:length(self.stats)
            temp = getfield(self.stats{k}, stat);
            fp(ax, self.stats{k}.wall_time+offset, temp, '-');
            hold(ax, 'on');
            fp(ax, self.stats{k}.wall_time(end)+offset, temp(end), 'o');
            offset = offset +self.stats{k}.wall_time(end);
        end
        grid(ax, 'on');
        xlabel(ax, 'wall time (s)')
        ylabel(ax, stat);
    end


    function plot_pop(self, varargin)
    % Plots the evolution of the populations under the current
    % control sequence.

        % what should we plot?
        if self.system.liouville
            plot_func = @prob_stateop;
        else
            plot_func = @prob_ket;
        end

        self.plot_X(plot_func, varargin{:});

        % NOTE: due to the horrible scoping rules of MATLAB, we use small x
        % in these functions as not to nuke the capital X in the parent function.

        function p = prob_stateop(x)
        % Returns the diagonal of a vectorized state operator.
            x = partial_trace(inv_vec(x), self.system.dimSE, 2);
            p = real(diag(x));
        end

        function p = prob_ket(x)
        % Returns the absolute values squared of ket vector elements.
            p = real(x .* conj(x));
        end
    end

    function plot_eig(self, varargin)

        % what should we plot?
        if self.system.liouville
            plot_func = @eig_stateop;
        else
            error('System closed, eigenvalues are constant.')
        end

        self.plot_X(plot_func, varargin{:});

        function p = eig_stateop(x)
        % Returns the real, nonnegative eigenvalues of the state operator.
            x = partial_trace(inv_vec(x), self.system.dimSE, 2);
            x = 0.5 * (x + x'); % eliminate numerical errors
            p = eig(x);
        end
    end

    function plot_X(self, plot_func, dt, k, ax, full_plot)
    % Plots the evolution of the initial system as a function of
    % time under the current control sequence.
    % dt, if given, is the timestep.
    % k is an optional ensemble index.
    % ax are the set of axes the plot goes into.
    % TODO for now it only handles kets and state ops

        if nargin < 6
            full_plot = true;
        if nargin < 5
            ax = gca();
        if nargin < 4
            k = 1; % TODO should match the choice in X()
        if nargin < 3
            dt = []; % one plot point per bin
        end
        end
        end
        end

        n_timeslots = self.seq.n_timeslots();
        n_ensemble  = self.system.n_ensemble();

        if full_plot
            % things that don't change and aren't deleted by cla
            set_plotstyle(ax);
            temp = self.system.description;
            if ~isempty(k)
                temp = sprintf('%s (%d)', temp, k);
            end
            title(ax, temp);
            xlabel(ax, 'time');
            ylabel(ax, 'probability');
            set(ax, 'NextPlot','replacechildren'); % so plot() won't reset these properties
        else
            %cla(ax);
        end

        if isempty(dt)
            % one plot point per timeslot
            for j = 0:n_timeslots
                res(j+1, :) = plot_func(self.X(j, k));
            end
            t = [0; cumsum(self.seq.tau)];
        else
            % use the given dt for plotting
            t = [0];
            X = self.X(0, k);
            res(1, :) = plot_func(X);

            for j = 1:n_timeslots
                X_end = self.X(j, k); % make sure the cache is up-to-date
                G = self.cache.H{j, k};
                tau = self.seq.tau(j);
                
                n_steps = ceil(tau / dt); % at least one point per timeslot
                step = tau / n_steps;
                P = expm(step * G);
                for q = 1:n_steps
                    X = P * X;
                    res(end+1, :) = plot_func(X);
                end
                temp = t(end);
                t = [t, linspace(temp+step, temp+tau, n_steps)];
                X = X_end; % stability...
            end
        end
        plot(ax, t, res);
        axis(ax, [0, t(end), 0, 1]);
        legend(self.system.state_labels, 'interpreter', 'latex');
    end

    function shake(self, rel_change)
    % Makes a small random perturbation to the control sequence, can be used
    % to shake the system out of a local optimum.
    % Does not touch the tau values.
        self.seq.shake(rel_change);
        self.cache.invalidate(); % flush the entire cache
    end

    function split(self, bins, n)
    % Refines the sequence by splitting the given bins into n equal pieces.
    % If an empty vector of bin numbers is given, the entire sequence is refined.
        self.seq.split(bins, n);
        self.cache_init(); % rebuild the cache
    end
  end
end
