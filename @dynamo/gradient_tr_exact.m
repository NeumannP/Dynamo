function ret = gradient_tr_exact(self, t, k, c)
% Exact gradient of error_tr.

% Uses the eigendecomposition.
%
% Uses H{t}, U{t_c}, U{t_tau+1} and L{t+1}.

if c < 0
    % dP_t/dtau_{t} = H_t P_t = P_t H_t
    dQdtau = partial_trace(self.cache.L{t+1, k} * self.cache.H{t, k} * self.cache.U{t+1, k}, self.system.dimSE, 1);
    ret = -self.seq.tau_deriv(t) * trace_matmul(self.cache.VUdagger, dQdtau);
else
    dPdu = dPdu_exact(self.cache.H_v{t, k}, self.cache.H_eig_factor{t, k}, self.system.B{k, c});
    dQdu = partial_trace(self.cache.L{t+1, k} * dPdu * self.cache.U{t, k}, self.system.dimSE, 1);
    ret = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * trace_matmul(self.cache.VUdagger, dQdu);
end
