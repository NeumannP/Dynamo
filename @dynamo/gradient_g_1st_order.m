function ret = gradient_g_1st_order(self, t, k, c)
% Gradient of the auxiliary function g by first order approximation.

% dP_k/du_c \approx (B_c * dt_k) * P_k
% Exact if G_k commutes with B_c.
%
% Uses H{t}, U{t+1} and L{t+1}.

if c < 0
    % tau control, exact
    ret = self.seq.tau_deriv(t) * trace_matmul(self.cache.L{t+1, k}, self.cache.H{t, k} * self.cache.U{t+1, k});
else
    % other control, approximate
    % TODO: this test is _really_ expensive, around 20% of total running time.
    %self.gradient_test(t, k, c);
    ret = self.seq.tau(t) * self.seq.fields_deriv(t, c) * trace_matmul(self.cache.L{t+1, k}, self.system.B{k, c} * self.cache.U{t+1, k});
end

ret = -ret * self.cache.VUdagger;
