classdef dynamo < matlab.mixin.Copyable
% Copyable handle class for DYNAMO optimizer objects.
%
% Contains the optimization task, system description, control
% sequence, various options and statistics. This class glues together
% the functionality of the cache, qsystem and control classes.
%
% Governing equation: \dot(X)(t) = (A +\sum_k u_k(t) B_k) X(t) = G(t) X(t)
    
% Shai Machnes   2010-2011
% Ville Bergholm 2011-2017


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
        ret = '1.4.0+';
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
    function self = dynamo(task, initial, final, A, B, weight, n_controls)
    % Constructor

        if nargin < 6
            % ensemble weights
            weight = 1;
        end
        if nargin < 7
            if iscell(B)
                if iscell(B{1})
                    % TEST TODO: allow nested cell vectors for the
                    % Python interface (which cannot do cell arrays).
                    ne = length(B);
                    nc = length(B{1});
                    BB = cell(ne, nc);
                    for k=1:ne
                        if length(B{k}) ~= nc
                            error('Each ensemble member must have the same number of controls.')
                        end
                        for c=1:nc
                            BB{k, c} = B{k}{c};
                        end
                    end
                    B = BB;
                end
                n_controls = size(B, 2);
            else
                error('If B is not a cell array, you must input the number of controls separately.')
            end
        end
        task = lower(task);

        %% Some basic data provenance

        config.version = dynamo.version();
        % Local time. UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
        config.date = datestr(now(), 31);
        % HACK: UTC time using Java
        temp = datenum('1970', 'yyyy') +java.lang.System.currentTimeMillis / 86400e3;
        config.date_UTC = datestr(temp, 31);
        
        %% Configuration defaults

        config.task = task;
        config.nonprojective_error = false;  % switch for error_abs.m
        config.frobenius_like_error = true;  % error can be interpreted as |A-B|_F^2 / (2 * |A|_F^2)
        config.epsilon  = 2e-8;   % finite differencing: approximately sqrt(eps(double))
        config.UL_mixed = false;  % HACK: mixed states in a closed system

        %% Description of the physical system

        input_dim  = [size(initial, 1), size(final, 1)];
        input_rank = [size(initial, 2), size(final, 2)]; % check the validity of the inputs
        sys = qsystem(weight, input_dim, n_controls);

        %% Optimization task

        [system_str, rem] = strtok(task);
        [task_str, rem] = strtok(rem);
        [extra_str, rem] = strtok(rem);

        out = 'Target operation:';
        switch system_str
          case 'abstract'
            %% No transformations done on the A and B operators. 
            % the generator may be anything, hence error_full
            config.dP = 'fd';
            config.error_func = @error_full;
            
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


          case {'closed'}
            %% Closed system
            % the generator is always Hermitian and thus normal => use exact gradient
            config.dP = 'eig';
            config.error_func = @error_abs;

            switch task_str
              case 'state'
                % TEST more efficient Hilbert space implementation
                out = strcat(out, ' mixed state transfer');
                if all(input_rank == 1)
                    error('Either the initial or the final state should be a state operator.')
                end
                sys.hilbert_representation(to_op(initial), to_op(final), A, B, false);
                config.UL_mixed = true;
                config.nonprojective_error = true;
                if strcmp(extra_str, 'overlap')
                    % overlap error function
                    % NOTE simpler error function and gradient, but final state needs to be pure
                    out = strcat(out, ' (overlap)');
                    config.f_max = 1;
                else
                    % full distance error function
                    config.f_max = (sys.norm2 +norm2(sys.X_initial)) / 2;
                end

              case {'ket', 'gate'}
                task_is_ket = strcmp(task_str, 'ket');
                if task_is_ket
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
                switch extra_str
                  case 'phase'
                    out = strcat(out, ' (with global phase (NOTE: unphysical!))');
                    config.nonprojective_error = true;
                  case 'subspace'
                    if task_is_ket
                        error('Subspace option only available for gates.')
                    end
                    out = strcat(out, sprintf(' (on a %d-dimensional subspace)', round(sys.norm2)));
                  otherwise
                    out = strcat(out, ' (ignoring global phase)');
                end

              % system S + environment E
              case 'state_partial'
                error('unfinished, gradient is complicated')
                out = strcat(out, ' partial mixed state transfer (on S)');
                sys.hilbert_representation(to_op(initial), to_op(final), A, B, false);
                config.error_func = @error_full;
                config.UL_mixed = true;

              case 'gate_partial'
                out = strcat(out, ' partial unitary gate (on S)');
                if any(input_rank == 1)
                    error('Initial and final states should be unitary operators.')
                end
                sys.hilbert_representation(initial, final, A, B, true);
                config.f_max = sys.norm2;
                config.error_func = @error_tr;
                
              otherwise
                error('Unknown task.')
            end
            out = strcat(out, ' in a closed system.\n');

            
          case {'open'}
            %% Open system with a Markovian bath
            % The generator isn't usually normal, so we cannot use the exact gradient method
            config.dP = 'fd';
            config.error_func = @error_full;

            switch task_str
              case 'state'
                out = strcat(out, ' state transfer');
                sys.vec_representation(initial, final, A, B, true);
                if strcmp(extra_str, 'overlap')
                    % overlap error function
                    % NOTE simpler error function and gradient, but final state needs to be pure
                    out = strcat(out, ' (overlap)');
                    config.error_func = @error_abs;
                    config.f_max = 1;
                else
                    % full distance error function
                end
                
              case 'gate'
                out = strcat(out, ' quantum gate');
                if any(input_rank == 1)
                    error('Initial and final states should be unitary operators.')
                end
                sys.vec_representation(initial, final, A, B, false);

              % system S + environment E
              case 'state_partial'
                out = strcat(out, ' partial state transfer (on S)');
                sys.vec_representation(initial, final, A, B, true);
                
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

        switch extra_str
          case 'overlap'
            if input_rank(2) ~= 1 || abs(sys.norm2 -1) > 1e-6
                error('Overlap error functions are only available if the target state is a normalized ket.')
            end
            config.nonprojective_error = true;
            config.frobenius_like_error = false;  % cannot be interpreted as Frobenius errors
          otherwise
        end

        fprintf(out);
        switch sys.type
          case 'liouville'
            fprintf('Liouville space dimension: %d\n', sys.dim^2);
          case 'abstract'
            fprintf('Vector space dimension: %d\n', sys.dim);
          otherwise
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

        % UL_mixed: mixed states in a closed system
        self.cache = cache(self.seq.n_timeslots(), self.system.n_ensemble(), U_start, L_end, self.config.dP, self.config.UL_mixed);
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

    function mask = set_state_transform(self, t, k, P)
    % Replaces the propagator for time slice t with an arbitrary state transform,
    % e.g. a subspace projector.
        self.cache.set_P(t, k, P);
        self.seq.raw(t, 1:end-1) = pi/2;  % set the phantom controls to zero to make the sequence look nicer
        self.seq.transform();
        % do not update the controls so as to not overwrite P_xy
        % TODO better solution needed
        mask = self.full_mask(true);
        mask(t,:) = false;
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


    function [err_out, grad_out] = compute_error(self, control_mask, ensemble_idx)
    % Returns the error (and its gradient) at current control values.
    % This is where we sum over the system ensemble if necessary.

        % gradient requires a control mask
        if nargout == 2 && (nargin < 2 || isempty(control_mask))
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
            % _tr, _abs:
            self.cache.g_needed_now = true;
        end
        
        self.cache_refresh(); % this call does the heavy computing (expm etc.)

        err_out  = 0;
        grad_out = 0;
        n_ensemble = self.system.n_ensemble();
        if nargin < 3
            % weighted average over the ensemble
            pick = 1:n_ensemble;
            weight = self.system.weight;
        else
            % just one ensemble member
            pick = ensemble_idx;
            weight = ones(1, n_ensemble);
        end
        for k = pick
            %% (real) normalized error
            err = self.config.error_func(self, k) / self.system.norm2;
            fprintf('Error (%d): %g\n', k, err);

            % weighted ensemble average
            err_out = err_out +weight(k) * err;
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
                grad(z) = self.config.error_func(self, k, t, c);
            end
            % real, normalized, weighted gradient
            grad_out = grad_out +(weight(k) / self.system.norm2) * grad;
        end
    end

    function [err] = frobenius_error(self, varargin)
    % Returns the Frobenius norm error at current control values.
    % Same input options as with the dynamo.X() method.
    % NOTE that this function does not return a meaningful result for all optimization tasks, due to global phase.

        err = norm(self.X(varargin{:}) -self.system.X_final, 'fro');
    end


    function update_controls(self, raw, control_mask)
    % Updates selected raw controls.
    %
    %  raw: vector of raw, untransformed control values, control_mask: corresponding mask.
    %
    %  Updates control values for which control_mask is true.
    %  Makes the changed timeslots and stuff that depends on them stale.

        old = self.seq.get_raw();
         
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
            self.seq.set_raw(new);
            self.cache.mark_as_stale(changed_t_mask);
        end
    end


    function set_controls(self, fields, tau)
    % Sets the controls to the given values.
    % If the requested values are incompatible with the current
    % control transforms, fails with an error message.

        n_timeslots = self.seq.n_timeslots();
        n_controls = self.seq.n_controls();

        % reverse transform into raw control parameters
        if isempty(fields)
            % use old values
            fields = self.seq.raw(:,1:end-1);
        else
            if isscalar(fields)
                % set all fields to the same value
                fields = fields * ones(1, n_controls);
            end
            fields = self.seq.inv_transform(fields);
            if size(fields, 1) == 1
                % set all timeslots to the same values
                fields = ones(n_timeslots, 1) * fields;
            end
        end
        if nargin < 3 || isempty(tau)
            % use old values
            tau = self.seq.raw(:,end);
        else
            tau = self.seq.inv_transform_tau(tau);
        end

        self.seq.set_raw([fields, tau]);
        self.cache.invalidate(); % flush the entire cache
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

        [stat, tail] = strtok(stat);
        [pt, tail] = strtok(tail);
        if isempty(tail)
            marker = '-';
        else
            marker = tail;
        end
        % plot type
        switch strtok(pt)
          case 'semilog'
            fp = @semilogy;
          otherwise
            fp = @plot;
        end
        switch stat
          case 'error'
            fp = @(a,t,x,m) fp(a, t, abs(x), m);
          case 'control_integral'
          otherwise
            error('Unknown stat.')
        end

        offset = 0;
        for k=1:length(self.stats)
            temp = getfield(self.stats{k}, stat);
            fp(ax, self.stats{k}.wall_time+offset, temp, marker);
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
         switch self.system.type
           case 'liouville'
             plot_func = @prob_stateop;
           case 'abstract'
             plot_func = @(x) x;
           otherwise
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
        switch self.system.type
          case 'liouville'
            plot_func = @eig_stateop;
          otherwise
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

    function plot_X(self, ev_func, dt, k, ax, full_plot)
    % Evaluates ev_func on the system state as a function of
    % time under the current control sequence, plots the result.
    % dt, if given, is the timestep.
    % k is an optional ensemble index.
    % ax are the set of axes the plot goes into.
        if nargin < 6
            full_plot = true;
        if nargin < 5
            ax = gca();
        if nargin < 4
            k = 1; % TODO should match the choice in X()
        if nargin < 3
            dt = [];
        end
        end
        end
        end
        
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

        [res, t] = self.evaluate_X(ev_func, dt, k);
        plot(ax, t, res);
        %axis(ax, [0, t(end), 0, 1]);
        legend(self.system.state_labels, 'interpreter', 'latex');
    end

    function [res, t] = evaluate_X(self, ev_func, dt, k)
    % Evaluates ev_func on the system state as a function of
    % time under the current control sequence.
    % dt, if given, is the timestep.
    % k is an optional ensemble index.

        if nargin < 4
            k = 1; % TODO should match the choice in X()
        if nargin < 3
            dt = []; % one point per bin
        end
        end
        n_timeslots = self.seq.n_timeslots();
        n_ensemble  = self.system.n_ensemble();
        res = [];
        if isempty(dt)
            % one point per timeslot
            for j = 0:n_timeslots
                res(j+1, :) = ev_func(self.X(j, k));
            end
            t = [0; cumsum(self.seq.tau)];
        else
            % use the given dt
            t = [0];
            X = self.X(0, k);
            res(1, :) = ev_func(X);

            for j = 1:n_timeslots
                X_end = self.X(j, k); % make sure the cache is up-to-date
                G = self.cache.H{j, k};
                tau = self.seq.tau(j);
                
                n_steps = ceil(tau / dt); % at least one point per timeslot
                step = tau / n_steps;
                P = expm(step * G);
                for q = 1:n_steps
                    if self.config.UL_mixed
                        X = P * X * P';
                    else
                        X = P * X;
                    end
                    res(end+1, :) = ev_func(X);
                end
                temp = t(end);
                t = [t, linspace(temp+step, temp+tau, n_steps)];
                X = X_end; % stability...
            end
        end
    end

    function shake(self, rel_change, abs_change, t_dependent)
    % Makes a small random perturbation to the control sequence, can be used
    % to shake the system out of a local optimum.
    % Does not touch the tau values.
        if nargin < 4
            t_dependent = false;
        end
        if nargin < 3
            abs_change = 0;
        end
        if nargin < 2
            rel_change = 0.1;
        end
        self.seq.shake(rel_change, abs_change, t_dependent);
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
