function ret = error_tr(self, k, t, c)
% Error function for closed S+E.
%
%  err  = error_tr(self, k)
%    E_tr(A,B) in ensemble k.
%
%  grad = error_tr(self, k, t, c)
%    Gradient of E_tr(A,B) in ensemble k, with respect to
%    control field c, at time slice t.


if nargin == 2
    % compute the error
    % trace norm
    [U, S, V] = svd(self.cache.g{k});
    ret = self.config.f_max -trace(S);

    % NOTE assumes S is not singular
    self.cache.VUdagger = V * U';
    self.cache.E = ret;
else
    % compute the gradient
    % Uses E and L{t+1}, and additionally whatever gradient_dPdu uses.

    if c > 0 && self.config.dP(1) == 'f'
        % We may compute the finite diff approximation at three locations: P, g, E
        % Here we handle the last two cases.
        % Note that finite diff makes no sense for a partial derivative
        % wrt. tau, since those can be computed analytically.

        [P, epsilon] = self.finite_diff_P(t, k, c);
        g = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 1);
        %dgdu = (g -self.cache.g{k}) / epsilon;
        %ret = -trace_matmul(self.cache.VUdagger, dgdu);
        %return
        E = self.config.f_max -sum(svd(g));
        ret = (E -self.cache.E) / epsilon;
        return
    end

    [dPdu, scale, ut] = self.gradient_dPdu(t, k, c);

    % auxiliary function g: trace over S
    dgdu = partial_trace(self.cache.L{t+1, k} * dPdu * self.cache.U{ut, k}, self.system.dimSE, 1);
    ret = -scale * trace_matmul(self.cache.VUdagger, dgdu);
    ret = real(ret);
end
