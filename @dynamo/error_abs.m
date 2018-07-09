function ret = error_abs(self, k, t, c)
% Error function and its gradient for closed systems.
%
%  Special case (dim E = 1) of error_tr().
%  Separately implemented for performance reasons.
%
%  err  = error_abs(self, k)
%    E_abs(A,B) in ensemble k.
%
%  grad = error_abs(self, k, t, c)
%    Gradient of E_abs(A,B) in ensemble k, with respect to
%    control field c, at time slice t.

% g(A,B) = \tr_S(A'*B) is enough to compute the error function
% D(A,B) whenever |B| is constant. This holds in every closed
% system task except "closed state_partial".


if nargin == 2
    % compute the error
    g = self.cache.g{k};

    if self.config.nonprojective_error
        ret = self.config.f_max -real(g);
        self.cache.VUdagger = 1; % HACK, no SVD involved
    else
        temp = abs(g);
        ret = self.config.f_max -temp;
        if temp == 0
            self.cache.VUdagger = 0; % |g| not differentiable at this point
        else
            self.cache.VUdagger = conj(g) / temp; % == V * U' for [U, S, V] = svd(g)
        end
    end
    self.cache.E = ret;
else
    % compute the gradient
    % Uses g, E and L{t+1}, and additionally whatever gradient_dPdu uses.

    if c > 0 && self.config.dP(1) == 'f'
        % We may compute the finite diff approximation at three locations: P, B_S, E
        % Here we handle the last two cases.
        % Note that finite diff makes no sense for a partial derivative
        % wrt. tau, since those can be computed analytically.

        [P, epsilon] = self.finite_diff_P(t, k, c);
        g = trace_matmul(self.cache.L{t+1, k}, P * self.cache.U{t, k});
        %ret = -(self.cache.VUdagger / epsilon) * (g -self.cache.g{k});
        %return
        if self.config.nonprojective_error
            E = self.config.f_max -real(g);
        else
            E = self.config.f_max -abs(g);
        end
        ret = (E -self.cache.E) / epsilon;
        return
    end

    [dPdu, scale, ut] = self.gradient_dPdu(t, k, c);

    if self.config.UL_mixed
        % Special case: mixed states in a closed system, 'eig' gradient, nonprojective error.
        % Uses P{t_c}, U{t+1}
        scale = scale * 2;
        ut = t+1;
        if c > 0
            dPdu = dPdu * self.cache.P{t,k}';
        end
    end

    ret = -self.cache.VUdagger * scale * trace_matmul(self.cache.L{t+1, k}, dPdu * self.cache.U{ut, k});
    ret = real(ret);
end
