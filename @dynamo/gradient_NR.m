function [J] = gradient_NR(self, t, k, c, Q)
% Exact Jacobian of the Newton-Raphson quality function L
% with respect to the controls specified in control_mask. Uses eigendecomposition.
% Returns a matrix.
% FIXME not functional

% compute the eigenvalue factors for the target function
[v, zeta] = eig_factors(Q, true);

% iterate over the indices of the true elements of control_mask
J = zeros(prod(size(self.system.X_final)), nnz(control_mask));

    if c == tau_c
        % dP_t/dtau_{t} = H_t P_t = P_t H_t
        temp = self.seq.tau_deriv(t) * self.cache.L{t+1} * self.cache.H{t} * self.cache.U{t+1};
    else
        % Compute the derivative dP_t/du_{tc} using the eigendecomposition of H(t)
        % For efficiency -tau(t) has been factored out and applied later.
        temp = self.cache.H_v{t}' * self.system.B{c} * self.cache.H_v{t};
        dPdu_eigenbasis = temp .* self.cache.H_eig_factor{t}; % note the elementwise multiplication .* here
        dPdu = self.cache.H_v{t} * dPdu_eigenbasis * self.cache.H_v{t}'; % to computational basis

        temp = self.seq.tau(t) * self.seq.fields_deriv(t, c) * self.cache.L{t+1} * dPdu * self.cache.U{t};
    end
    J(:, z) = vec(project_to_su(v * ((v' * temp * v) ./ zeta) * v'));
end
