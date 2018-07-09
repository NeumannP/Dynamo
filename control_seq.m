classdef control_seq < matlab.mixin.Copyable
% Copyable handle class for control sequences.

% Ville Bergholm 2011-2016

  properties
      tau_par        % parametrization of tau, size == [n_timeslots, 2]
      control_type   % char array denoting the control type, size == [n_controls]
      control_par    % cell array of control parameter structs, size == [n_controls]

      raw            % untransformed control fields, last column is for tau: size == [n_timeslots, n_controls + 1]
      % the data below are computed from raw using the control parametrization
      tau            % timeslice duration, size == [n_timeslots]
      tau_deriv      % d tau / d tau_raw
      fields         % transformed control fields, size == [n_timeslots, n_controls]
      fields_deriv   % d fields / d fields_raw
  end

  methods
    function self = control_seq(n_timeslots, n_controls, tau_par, control_type, control_par)
    %  tau_par: [T_min, T_delta], either a single entry with the totals or one such entry for each time slot.

        %% Check control types
        if nargin < 4
            control_type = char(zeros(1, n_controls) + '.');
        elseif length(control_type) ~= n_controls
            error('Length of the control_type vector does not match the number of controls.')
        end

        fprintf('Timeslots: %d\nControls per slot: %d + tau\n\n', n_timeslots, n_controls);
        

        %% Check control parameters
        % For now control_par doesn't have separate entries for each timeslot
        if nargin < 5
            control_par = {};
        end

        %% Check tau parameters
        % Is tau_par a total time or a vector of delta_t:s?
        len_params = size(tau_par, 1);
        if len_params == 1
            % divide the parameters uniformly into timeslots
            tau_par = ones(n_timeslots, 1) * (tau_par / n_timeslots);
        elseif len_params ~= n_timeslots
            error('Length of tau_par does not match the number of control timeslots given.')
        end

        self.tau_par = tau_par;
        self.control_type = control_type;
        self.control_par = control_par;

        % init raw params to some reasonable values, but do not transform them yet
        self.raw = [zeros(n_timeslots, n_controls), acos(0) * ones(n_timeslots, 1)]; % tau: halfway
    end

    
    function ret = n_timeslots(self)
    % Returns the number of time slots.
        ret = size(self.tau_par, 1);
    end

    
    function ret = n_controls(self)
    % Returns the number of control fields (not including tau).
        ret = length(self.control_type);
    end
    
    
    function ret = get_raw(self, control_mask)
    % Returns the raw controls corresponding to the mask given, or all
    % of them if no mask is given.

        if nargin == 1
            ret = self.raw;
        else
            ret = self.raw(control_mask);
        end
    end

    
    function set_raw(self, raw)
    % Transform and set the controls using a diagonal transformation function.
    %
    % raw: raw, untransformed control values, size(raw) == [n_timeslots, n_controls + 1].
    % Given raw controls r_k, computes the transformed controls c_k(r_k) and their derivatives.
    % Uses the (fixed) control parameters.

        sss = size(raw);
        
        % Check the number of control fields. The last column are the tau controls.
        if sss(1) ~= self.n_timeslots()
            error('Given controls have the wrong number of timeslots.')
        end
        if sss(2) ~= self.n_controls() + 1
            error('Given controls have the wrong number of control fields.')
        end
        self.raw = raw;
        self.transform();
    end

    function transform(self)
    % Transforms the array of raw controls, and use it to initialize
    % taus and control fields and their derivatives.

        % Tau control is the last one.
        t_raw = self.raw(:, end);

        % Returned tau values should always be positive, and not too large
        % if we're using the 1st order gradient approximation.
        % Stretchy bins with min and max duration:
        % t_raw = 0    <=> min duration
        % t_raw = pi/2 <=> average duration
        self.tau       = self.tau_par(:,1) +0.5 * self.tau_par(:,2) .* (1-cos(t_raw));
        self.tau_deriv = 0.5 * self.tau_par(:,2) .* sin(t_raw);

        sss = size(self.raw) -[0, 1];
        self.fields = zeros(sss);
        self.fields_deriv = zeros(sss);

        for k=1:self.n_controls()
            temp = self.raw(:, k);
            switch self.control_type(k)
              case '.'  % no transformation
                self.fields(:, k) = temp;
                self.fields_deriv(:, k) = 1;
        
              case 'p'  % strictly nonnegative, u_k = r_k^2
                self.fields(:, k) = temp.^2;
                self.fields_deriv(:, k) = 2 * temp;
                
              case 'm'  % minimum and delta, u_k = min + delta * 0.5 * (1 - cos(r_k))
                par = self.control_par{k};
                self.fields(:, k) = par(1) +par(2) * 0.5 * (1 -cos(temp));
                self.fields_deriv(:, k) = par(2) * 0.5 * sin(temp);
      
              otherwise
                error('Unknown control type.')
            end
        end
        % TODO u_x = A*cos(phi), u_y = A*sin(phi) (not diagonal, but within the same bin, i.e. block diagonal)
        % fields_deriv(slot, c_j, raw_k) = d c_{sj} / d raw_{sk}
    end

    function ret = inv_transform(self, fields)
    % Given an array of control field values (columns corresponding
    % to different fields, not including tau), returns the corresponding raw controls.

        n_controls = self.n_controls();
        if size(fields, 2) ~= n_controls
            error('Wrong number of control fields.')
        end

        ret = zeros(size(fields));
        for k=1:n_controls
            switch self.control_type(k)
              case '.'  % no transformation
                ret(:, k) = fields(:, k);
        
              case 'p'  % strictly nonnegative, u_k = r_k^2
                if any(fields(:, k) < 0)
                    error('Field %d not nonnegative.', k)
                end
                ret(:, k) = sqrt(fields(:, k));
                
              case 'm'  % minimum and delta, u_k = min + delta * 0.5 * (1 - cos(r_k))
                par = self.control_par{k};
                ret(:, k) = acos(1 -(fields(:, k) -par(1)) * (2 / par(2)));
                if any(imag(ret(:, k)))
                    warning('Field %d not within the parameter limits.', k)
                    crap = imag(ret(:, k))
                    ret(:, k) = real(ret(:, k));
                end
              otherwise
                error('Unknown control type.')
            end
        end
    end

    
    function ret = inv_transform_tau(self, tau)
    % Given a column vector of tau values returns the corresponding raw tau controls.
    
        if size(tau, 1) ~= self.n_timeslots()
            error('Wrong number of time slots.')
        end

        ret = acos(1 -(tau -self.tau_par(:,1)) .* (2 ./ self.tau_par(:,2)));
        if any(imag(ret))
            warning('Some taus not within the parameter limits.')
            crap = imag(ret)
            ret = real(ret);
        end
    end
    
    
    function split(self, bins, n)
    % Refines the sequence by splitting the given bins into n equal pieces.
    % If an empty vector of bin numbers is given, the entire sequence is refined.

        % see how many timeslots we have after the split
        n_timeslots_old = self.n_timeslots();
        if isempty(bins)
            bins = 1:n_timeslots_old;
        else
            bins = unique(bins); % also sorts
        end
        n_timeslots_new = n_timeslots_old + (n-1) * length(bins);

        raw = zeros(n_timeslots_new, size(self.raw, 2));
        tau_par = zeros(n_timeslots_new, 2);
        % source and destination indices (first unprocessed slots)
        si = 1;
        di = 1;
        for k=1:length(bins)
            b = bins(k);
            % a run of unchanged slots
            run = b -si -1; % run length-1
            tau_par(di:di+run, :) = self.tau_par(si:si+run, :); 
            raw(di:di+run, :) = self.raw(si:si+run, :); 
            di = di+run+1;
            % a split slot
            % tau_par slots are actually split, raw slots are just duplicated
            % this way we don't have to care about inverse transforms here.
            % TODO should the split bins be able to strech as large as the parent?
            tau_par(di:di+n-1, :) = ones(n, 1) * self.tau_par(b, :) / n;
            raw(di:di+n-1, :) = ones(n, 1) * self.raw(b, :);
            si = b + 1;
            di = di + n;
        end
        % final run of unchanged slots
        tau_par(di:end, :) = self.tau_par(si:end, :); 
        raw(di:end, :) = self.raw(si:end, :); 

        % transform the new controls
        self.tau_par = tau_par;
        self.raw = raw;
        self.transform();
    end


    function clip(self, limit)
    % Clips the controls according to the limits in the limit vector.

        n_timeslots = self.n_timeslots();
        n_controls = self.n_controls();
        t_shape = [n_timeslots, 1];
        limit = ones(t_shape) * limit;

        f = self.raw(:, 1:n_controls);
        mask = f > limit;
        f(mask) = limit(mask);
        mask = f < -limit;
        f(mask) = -limit(mask);
        self.raw(:, 1:n_controls) = f;
        self.transform();
    end

    function shake(self, rel_change, abs_change, t_dependent)
    % Makes a small random perturbation to the control sequence, can be used
    % to shake the system out of a local optimum.
    % Does not touch the tau values.

        n_timeslots = self.n_timeslots();
        n_controls = self.n_controls();
        % shape vectors
        f_shape = [n_timeslots, n_controls];
        t_shape = [n_timeslots, 1];
        c_shape = [1, n_controls];

        f = self.raw(:, 1:n_controls);
        if ~t_dependent
            % perturb each control field randomly
            p_abs = ones(t_shape) * (abs_change .* randn(c_shape));
            p_rel = ones(t_shape) * (1 +rel_change .* randn(c_shape));
            f = f .* p_rel +p_abs;
        else
            % add some time-dependent noise on top of the raw controls
            f = f .* (1 +rel_change * randn(f_shape)) +abs_change * randn(f_shape);
        end
        self.raw(:, 1:n_controls) = f;
        self.transform();
    end


    function ret = fluence(self, M)
    % Computes the total fluence of the control fields.
    %
    % The fluence is defined as F^2 = \int_0^T |H_c(t)|^2 dt,    H_c(t) = \sum_r H_r f_r(t).
    % If the control Hamiltonians are orthonormal, M_ij = (H_i, H_j) = \delta_{ij},
    % this gives F^2 = \sum_k \int_0^T |f_k(t)|^2 dt.
    %
    % For piecewise constant controls this reduces to F^2 = \sum_{ri} |a_{ri}|^2 \Delta t_{i}
    % F^2 has units of energy^2 * time.
    % TODO relation to signal energy? E = F^2 / \hbar?

        ret = 0;
        n_timeslots = self.n_timeslots();
        for k = 1:n_timeslots
            c = self.fields(k, :); % all controls for this timeslot
            ret = ret +c * M * c.' * self.tau(k);
        end
    end

    
    
    function ret = integral(self)
    % Computes the time integral of the control fields.

        ret = transpose(self.tau) * self.fields;
    end
    


    function ax = plot(self, labels, ax, full)
    % Plots a control sequence using superposed bars.

        if nargin < 3
            % If no axes are given, use the current ones.
            ax = gca();
        end
        if nargin < 4
            full = true;
        end
        
        nc = self.n_controls();
        
        if full
            % things that don't change and aren't deleted by cla
            set_plotstyle(ax);
            colormap(jet)
            set(ax, 'CLim', [0 nc])
            title(ax, 'Control Sequence')
            xlabel(ax, 'time')
            ylabel(ax, 'control amplitudes')
        else
            cla(ax);
        end
        
        % start times for pulses
        t = [0; cumsum(self.tau)];
        c = self.fields;

        % set new axis limits
        axis(ax, [0, t(end), min(min(c(:)), 0), max(max(c(:))+1e-3, 0)]);
        hold(ax, 'on');
        for j=1:length(t)-1
            x = [t(j), t(j+1), t(j+1), t(j)];
    
            % use painter's algorithm
            [dummy, I] = sort(abs(c(j, :)), 2, 'descend'); % dummy instead of ~ for Octave
            for k=1:nc
                y = c(j, I(k));
                p(k) = patch(x, [0, 0, y, y], 1);
                set(p(k), 'FaceColor','flat',  'CData',I(k), 'CDataMapping','scaled')
            end
            if j == 1
                p_colors(I) = p; % HACK, to get the legend right
            end
        end
        hold(ax, 'off')
        legend(ax, p_colors, labels); % cla deletes this, so it needs to be redrawn

        %stairs(t, c)
        %c = [c; zeros(1,nc)]; % final time step is a dummy
        %for j=1:k
        %    bar(t, c(:,j), 1, 'histc')
        %end
    end
  end
end

