classdef cache < matlab.mixin.Copyable
% Copyable handle class for doing the heavy computing and storing the results.

% Shai Machnes   2010-2011
% Ville Bergholm 2011-2013
    

  properties (SetAccess = private)
      %% cell arrays: Q{timeslice, ensemble_index}
      H  % Generator for a time slice. Not just the Hamiltonian, see below.
      P  % Propagator for a time slice, computed by calcPfromHfunc. With a constant H, P == expm(dt * H). 
      U  % Forward propagators. U{k+1} = P{k} * U{k}
         % U{k} is the system at t = sum(tau(1:(k-1))) = t_{k-1}
      L  % Backward propagators. L{k-1} = L{k} * P{k-1};
         % L{k} is the adjoint system at t = sum(tau(1:(k-1))) = t_{k-1}

      H_v           % eigendecomposition data for dt*H, updated when P is updated  
      H_eig_factor  % likewise
      
      %% cell array: g{ensemble_index}
      g  % trace_q(L{k} * U{k}), where q is S or E depending on the error function used.
  end

  properties (Access = public)
      H_needed_now  % flags, Q[timeslice]
      P_needed_now
      U_needed_now
      L_needed_now
      g_needed_now  % scalar
      
      E        = NaN  % gradient_*_finite_diff: cached error
      VUdagger = NaN  % cached SVD data for error_tr, error_abs (also used by error_real) 
      int  = []   % TEST, integrator
  end

  properties (Access = private)
      H_is_stale  % flags, Q[timeslice]
      P_is_stale
      U_is_stale
      L_is_stale
      g_is_stale  % scalar
      
      calcPfromHfunc
      UL_mixed  % flag: use mixed states in Hilbert space representation?
  end

  methods
      function self = cache(n_timeslots, n_ensemble, U_start, L_end, use_eig, UL_hack)
      % Set up caching (once we know the number of time slices and U and L endpoints).

          s_time     = [1, n_timeslots];
          s_ensemble = [1, n_ensemble];
          s_full     = [n_timeslots, n_ensemble];
          
          if use_eig
              % Store the eigendecomposition data as well
              self.H_v          = cell(s_full);
              self.H_eig_factor = cell(s_full);
              self.calcPfromHfunc = @calcP_expm_exact_gradient;
          else
              self.calcPfromHfunc = @calcP_expm;
          end
          self.UL_mixed = UL_hack;

          self.H = cell(s_full);
          self.P = cell(s_full);
          self.U = cell(s_full + [1, 0]); % one more timeslot
          self.L = cell(s_full + [1, 0]); % one more timeslot

          for k=1:n_ensemble
              % U: X_initial propagated forward up to a time instant.
              self.U{1, k} = U_start;
              % L: X_final' propagated backward, or in open-system tasks, a pure propagator (even though this requires more memory).
              self.L{end, k} = L_end;
          end

          self.g = cell(s_ensemble);
          
          % Keep track of which objects need re-computation.
          % U{1} and L{end} are never stale and never recomputed.
          self.H_is_stale = true(s_time);
          self.P_is_stale = true(s_time);
          self.U_is_stale = [false, true(s_time)]; % Updates for H via 'control_update' get propagated automatically
          self.L_is_stale = [true(s_time), false];
          self.g_is_stale = true;
          
          % Here we indicate which objects we need to have up-to-date
          % The 'needed_now' function looks at everything we need to recompute now, 
          % compared to what in principle is_stale, and (in principle) executes the
          % optimal set of operations so that for everything which was marked
          % 'needed_now' is up-to-date

          % Example: We modified H{3}, so in theory we need to recompute U{4:end}.
          % But for our immediate needs we only want U{7} and L{7},
          % so we mark U{4:end} as "is_stale", but only 4:7 as "needed_now".

          self.H_needed_now = false(s_time);   
          self.P_needed_now = false(s_time);   
          self.U_needed_now = [false, false(s_time)];
          self.L_needed_now = [false(s_time), false]; 
          self.g_needed_now = false;
      end
      
      
      function invalidate(self)
      % Invalidates the entire cache.

          self.H_is_stale(:) = true;
          self.P_is_stale(:) = true;
          % the first U and last L never become stale
          self.U_is_stale(2:end) = true; 
          self.L_is_stale(1:(end-1)) = true;
          self.g_is_stale = true;
      end


      function mark_as_stale(self, changed_t_mask)
      % Marks the selected timeslots as stale.

          self.H_is_stale(changed_t_mask) = true;
          self.P_is_stale(changed_t_mask) = true;

          % Propagate the H_is_stale to the U and Ls.
          self.U_is_stale( (find(self.H_is_stale, 1, 'first')+1):end) = true;
          self.L_is_stale(1:find(self.H_is_stale, 1, 'last'))         = true;
          self.g_is_stale = true;
      end


      function set_P(self, t, k, P)
      % Sets P{t, k} to the given value.
          self.P{t, k} = P;
          self.H{t, k} = NaN;
          % mark U, L and g as stale, P and H as not stale.
          self.mark_as_stale(t);
          self.P_is_stale(t) = false;
          self.H_is_stale(t) = false;
      end


      function refresh(self, sys, tau, fields)
      % This function does most of the heavy computing.
      % It should be called _after_ setting _all_ the required *_needed_now fields
      % to minimize unnecessary computational work.
      %
      % The function looks at everything we need to recompute now,
      % compared to what in principle is_stale, and (in principle) executes the
      % optimal set of operations so that for everything which was marked
      % 'needed_now' is up-to-date.
      %
      % It then updates the _is_stale flags to indicate what has been actually
      % computed, and clears the _needed_now flags.
      %
      % The inter-dependence of H and U/L updates is taken care of in mark_as_stale()

          n_timeslots = size(self.H, 1);
          n_ensemble  = size(self.H, 2);
          n_controls  = sys.n_controls();
          
          % computing g may require additional U and L elements, so check that first
          g_recompute_now = self.g_needed_now && self.g_is_stale;
          if g_recompute_now
              % g can be computed using any slice k \in [1, n+1]: g = trace_S(L_k * U_k).
              % Try to figure out which k requires least additional computation.
              g_ind = self.g_setup_recalc();
          end
          
          U_recompute_now = self.U_needed_now & self.U_is_stale;
          L_recompute_now = self.L_needed_now & self.L_is_stale;
          % To recompute U, you need to start at a cell that is fully recomputed.
          % U{1} and L{n_timeslots+1} are by definition never stale and never recomputed.
          for t=n_timeslots:(-1):2
              if U_recompute_now(t+1) && self.U_is_stale(t)
                  U_recompute_now(t) = true;
              end
          end
          for t=2:n_timeslots
              if L_recompute_now(t-1) && self.L_is_stale(t)
                  L_recompute_now(t) = true;
              end
          end

          % Now that we know which Us & Ls we need to recompute, we can figure out which Ps and Hs must be up-to-date
          P_recompute_now = (U_recompute_now(2:end) | L_recompute_now(1:(end-1)) | self.P_needed_now) & self.P_is_stale;
          H_recompute_now = (P_recompute_now | self.H_needed_now) & self.H_is_stale;

          % timeslot indices to recompute
          h_idx = find(H_recompute_now);
          p_idx = find(P_recompute_now);
          u_idx = find(U_recompute_now);
          el_idx = fliplr(find(L_recompute_now));

      % loop over the ensemble of systems
      for k=1:n_ensemble

          if isequal(self.calcPfromHfunc, @calcP_int)
              % set up stuff in the integrator

              self.int.import_better(tau, fields, false, true);

              % integrate props.
              
              % TODO FIXME crosstalk makes all gradients except the
              % finite_diff ones invalid?, and even the finite_diff
              % ones have to be computed using an integrator?

              % FIXME with crosstalk currently taus must not change
              % (because a changing tau would invalidate every
              % propagator following it too)
          end
          % Compute the Hamiltonians
          for t=h_idx
              H = sys.A{k};
              for c = 1:n_controls
                  u = fields(t, c);
                  H = H +u * sys.B{k, c};
              end
              self.H{t, k} = H;
          end

          % Compute the exp(H) and any other per-H computation which may be needed for the gradient function
          for t=p_idx
              self.P{t, k} = self.calcPfromHfunc(self, t, k, tau(t)); % Compute the Ps - a single piece of propagator
              % NOTE: calcPfromHfunc may also compute other values which will be needed for gradient calculations.
              %       These should be stored in cache. Their up-to-date-ness is identical to that of P.
          end

          % Compute the Us - forward propagation (we never recompute U{1})
          % Compute the Ls - adjoint system propagation
          if self.UL_mixed
              % mixed states, unitary evolution: propagate from both sides
              for t=u_idx
                  self.U{t, k} = self.P{t-1, k} * self.U{t-1, k} * self.P{t-1, k}';
              end
              for t=el_idx
                  self.L{t, k} = self.P{t, k}' * self.L{t+1, k} * self.P{t, k};
              end
          else
              % propagate U from left, L from right
              for t=u_idx
                  self.U{t, k} = self.P{t-1, k} * self.U{t-1, k};
              end
              for t=el_idx
                  self.L{t, k} = self.L{t+1, k} * self.P{t, k};
              end
          end

          % and finally g := trace_q(L{g_ind} * U{g_ind})
          if g_recompute_now
              if self.g_needed_now == 2  % HACK ln37ae983e
                  % error_full, L:s are full propagators, partial trace over E gives X_S
                  temp = self.L{g_ind, k} * self.U{g_ind, k};
                  self.g{k} = partial_trace(temp, sys.dimSE, 2);
              
              elseif sys.dimSE(2) == 1
                  % error_abs, error_real, no environment E, full trace, g is a scalar
                  self.g{k} = trace_matmul(self.L{g_ind, k}, self.U{g_ind, k});
              
              else
                  % error_tr, partial trace over S, g is a matrix
                  temp = self.L{g_ind, k} * self.U{g_ind, k};
                  self.g{k} = partial_trace(temp, sys.dimSE, 1);
              end
              % TODO combine last 2 cases somehow?
          end

      end % loop over ensemble
              
          % Mark what has been actually computed
          temp = [1, n_timeslots];

          self.H_is_stale(H_recompute_now) = false;
          self.P_is_stale(P_recompute_now) = false;
          self.U_is_stale(U_recompute_now) = false;
          self.L_is_stale(L_recompute_now) = false;
          if g_recompute_now
              self.g_is_stale = false;
          end
          
          self.H_needed_now = false(temp);
          self.P_needed_now = false(temp);
          self.U_needed_now = [false, false(temp)];
          self.L_needed_now = [false(temp), false];
          self.g_needed_now = false;

          %ret = [sum(H_recompute_now), sum(P_recompute_now), sum(U_recompute_now), sum(L_recompute_now)]
      end


      function ret = calcP_expm_exact_gradient(self, t, k, dt)
      % Computes P{t, k} using the eigendecomposition, stores some
      % extra stuff for cheap exact gradient computation later on.

          dt_H = dt * self.H{t, k};
          %N = length(dt_H);

          %% Compute the eigenvalue factors and eigenvectors of dt*H

          [v, zeta, exp_d] = eig_factors(dt_H, true);
          self.H_v{t, k} = v;
          self.H_eig_factor{t, k} = zeta;

          %% And finally expm(dt*H) using the eigendecomposition

          ret = v * diag(exp_d) * v';
      end


      function ret = calcP_expm(self, t, k, dt)
      % Computes P{t, k} using expm.
          ret = expm(dt * self.H{t, k});
      end


      function ret = calcP_int(self, t, k, dt)
      % Computes P{t, k} by integrating a time-dependent Hamiltonian. 
      % Used to deal with control crosstalk in the RWA. Slow.

          % the integrator already has the fields and taus updated
          ret = self.int.int_bin_rwa(t);
      end

      function switch_to_int(self, int)
      % FIXME test, switch to crosstalk-aware integrator.
          self.int = int;
          self.calcPfromHfunc = @calcP_int;
          
          % force P, U, L recomputation
          self.invalidate();
      end
      

      function t = g_setup_recalc(self)
      % Returns the optimal time slice in which to compute g.

      % TODO: Optimally, we can search for the time which will require minimum calculations to compute the value.
      % However, this is tricky (we need to count the expensive H-->P computations (expm or similar) and maybe weigh-in the cheap
      % U and L updates (matrix multiplications).

      % But for now, we'll do a sub-optimal algorithm and take the left-most L that is up-to-date (and therefore update all L-s to the
      % right of it). If all the slots updated are in a single block, we'll be optimal

          flags = (~self.U_is_stale) + (~self.L_is_stale)*2;

          flag3 = find(flags==3,1,'first');
          if ~isempty(flag3)
              t = flag3;
          else
              p1 = find(flags>0);
              p2 = find((flags(p1(1:(end-1)))==1) & (flags(p1(2:end))==2));

              if isempty(p2)
                  error('There should be at least one');
              end
              cost = p1(p2+1)-p1(p2);
              %    [~,mincostpos] = min(cost);
              [dummy, mincostpos] = min(cost); % for Octave
              t = p1(p2(mincostpos));
          end
          self.U_needed_now(t) = true;
          self.L_needed_now(t) = true;
      end
  end
end
