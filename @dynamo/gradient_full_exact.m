function ret = gradient_full_exact(self, t, k, c)
% Gradient of error_full by auxiliary matrix method (exact).

% expm([G, B_c; 0, G]*t) = [P, dP/du_c; 0, P]
%
% Uses H{t}, U{t_c}, U{t_tau+1} and L{t+1}.

X_S = self.cache.g{k};
temp = X_S -self.system.X_final;

if c < 0
    % tau control
    ret = self.seq.tau_deriv(t) * inprod(temp, partial_trace(self.cache.L{t+1, k} * self.cache.H{t, k} * self.cache.U{t+1, k}, self.system.dimSE, 2));
else
    % other control
    d = self.system.dimSE(1);
    G = self.cache.H{t, k};
    mat = expm([G, self.system.B{k, c}; zeros(d), G] * self.seq.tau(t));
    dPdu = mat(1:d, d+1:end);
    ret = self.seq.fields_deriv(t, c) * inprod(temp, partial_trace(self.cache.L{t+1, k} * dPdu * self.cache.U{t, k}, self.system.dimSE, 2));
end
