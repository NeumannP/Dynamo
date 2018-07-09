function gradient_test(self, t, k, c)
% Checks whether the gradient commutator series has trouble converging
% at time slice t, ensemble member k, control c.

% The approximation is exact when this vanishes.
violation = abs(self.seq.tau(t)) * norm(self.cache.H{t, k});

if violation > self.opt.max_violation
    self.opt.max_violation = violation;
    
    if violation > 1    
        fprintf('warning: gradient approximation not valid at t = %d, k = %d, c = %d: violation = %f.\n', t, k, c, violation)
    end
end
