function [P, epsilon] = finite_diff_P(self, t, k, c)
% Computes the propagator P_{t, k} with the raw control c increased by epsilon.

% Uses H{t}.
    
epsilon = self.config.epsilon;


if isempty(self.cache.int)
    % FIXME integrator test
    
if c < 0
    tau_eps = self.seq.tau(t) +self.seq.tau_deriv(t) * epsilon;
    P = expm(tau_eps * self.cache.H{t, k});
else
    H_eps = self.cache.H{t, k} +(epsilon * self.seq.fields_deriv(t, c)) * self.system.B{k, c};
    P = expm(self.seq.tau(t) * H_eps);
end

else
    %cache.with_modified_control(t,c, f'(t,c)).integrate_bin(t, k);
    % f'(t,c) \approx f(t,c) +epsilon * self.seq.fields_deriv(t, c)
    % FIXME crosstalk!
    
    if c < 0
        % NOTE with crosstalk and varying taus, all the bins after the one whose
        % length we change will have different propagators because the
        % phases of the slowly rotating terms depend on (absolute) time! keep taus constant?
        error('do not do this')
    else
        ccc = mod(c-1, 2); % 0 or 1, x or y?
        temp = c -ccc;
        cn = 1 +(temp -1)/2; % carrier number

        % get the current x and y control pair
        a = self.seq.fields(t, temp:temp+1);
        % update them
        a(ccc+1) = a(ccc+1) +epsilon;
        
        amp = self.cache.int.amp(t, cn);
        phi = self.cache.int.phi(t, cn);
        
        % transfer to integrator
        self.cache.int.amp(t, cn) = sqrt(a(1)^2 +a(2)^2);
        self.cache.int.phi(t, cn) = atan2(a(2), a(1));
              
        P = self.cache.calcP_int(t, k);
        
        % restore originals
        self.cache.int.amp(t, cn) = amp;
        self.cache.int.phi(t, cn) = phi;
    end
end