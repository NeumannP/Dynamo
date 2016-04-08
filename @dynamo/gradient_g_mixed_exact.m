function ret = gradient_g_mixed_exact(self, t, k, c)
% Exact gradient of the auxiliary function g for mixed states in a
% closed system.

% Gradient of g with respect to the controls
% specified in control_mask. Uses the eigendecomposition.
%
% Uses H{t}, P{t}, U{t_c}, U{t_tau+1} and L{t+1}.

if c < 0
    % dP_t/dtau_{t} = H_t P_t = P_t H_t
    ret = 2 * self.seq.tau_deriv(t) * trace_matmul(self.cache.L{t+1, k}, self.cache.H{t, k} * self.cache.U{t+1, k} * self.cache.P{t, k}');
else
    dPdu = dPdu_exact(self.cache.H_v{t, k}, self.cache.H_eig_factor{t, k}, self.system.B{k, c});
    ret = 2 * self.seq.tau(t) * self.seq.fields_deriv(t, c) * ...
          trace_matmul(self.cache.L{t+1, k}, dPdu * self.cache.U{t, k} * self.cache.P{t, k}');
end

% HACK this is a single-use function (only used with error_real) so
% no VUdagger is needed, we can do the -1:s here directly
ret = -ret;
