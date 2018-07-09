function ret = error_full(self, k, t, c)
% Error function and its gradient for open (Markovian) systems.
%
%  err  = error_full(self, k)
%    E_full(A,B) in ensemble k.
%
%  grad = error_full(self, k, t, c)
%    Gradient of E_full(A,B) in ensemble k, with respect to
%    control field c, at time slice t.


if nargin == 2
    % compute the error
    temp = self.cache.g{k} -self.system.X_final;
    ret = 0.5 * norm2(temp);

    self.cache.E = ret;
else
    % compute the gradient
    % Uses g, E and L{t+1}, and additionally whatever gradient_dPdu uses.

    % \tr_E(B) -A
    delta = self.cache.g{k} -self.system.X_final;

    if c > 0 && self.config.dP(1) == 'f'
        % We may compute the finite diff approximation at three locations: P, B_S, E
        % Here we handle the last two cases.
        % Note that finite diff makes no sense for a partial derivative
        % wrt. tau, since those can be computed analytically.

        [P, epsilon] = self.finite_diff_P(t, k, c);
        B_S = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 2);
        %dB_Sdu = (B_S -self.cache.g{k}) / epsilon;
        %ret = inprod(delta, dB_Sdu);
        %return
        E = 0.5 * norm2(B_S -self.system.X_final);
        ret = (E -self.cache.E) / epsilon;
        return
    end

    [dPdu, scale, ut] = self.gradient_dPdu(t, k, c);
    ret = scale * inprod(delta, partial_trace(self.cache.L{t+1, k} * dPdu * self.cache.U{ut, k}, self.system.dimSE, 2));
    ret = real(ret);
end
