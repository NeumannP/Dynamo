function ret = gradient_full_1st_order(self, t, k, c)
% Gradient of error_full by first order approximation.

% dP_k/du_c \approx (B_c * dt_k) * P_k
% Exact if G_k commutes with B_c.
%
% Uses H{t}, U{t+1} and L{t+1}.

X_S = self.cache.g{k};
temp = X_S -self.system.X_final;
    
if c < 0
    % tau control, exact
    ret = self.seq.tau_deriv(t) * inprod(temp, partial_trace(self.cache.L{t+1, k} * self.cache.H{t, k} * self.cache.U{t+1, k}, self.system.dimSE, 2));
else
    % other control, approximate
    %self.gradient_test(t, k, c);
    ret = self.seq.tau(t) * self.seq.fields_deriv(t, c) * inprod(temp, partial_trace(self.cache.L{t+1, k} * self.system.B{k, c} * self.cache.U{t+1, k}, self.system.dimSE, 2));
end
