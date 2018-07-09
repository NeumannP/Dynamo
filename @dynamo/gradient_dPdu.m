function [dPdu, scale, ut] = gradient_dPdu(self, t, k, c)
% Partial derivative of P{t, k} with respect to u_c.

% Uses
% U{t}   fd, fd_dp, eig, aux, series_ss
% U{t+1} tau, series (implies P{t})
% P{t}   fd_dp, eig, series_ss, (implies H{t})
% H{t}   tau, fd, aux, series, series_ss

if c < 0
    % tau control:  dP/dtau = G P = P G  (exact)
    dPdu  = self.cache.H{t, k};
    scale = self.seq.tau_deriv(t);
    ut = t+1;  % missing P{t} for dPdu
else
    % normal control
    ut = t;

    switch self.config.dP
      case 'fd'
        % f'(x) = (f(x + eps) -f(x))/eps
        % Trivial and relatively slow, but a good reference point.
        [P, epsilon] = self.finite_diff_P(t, k, c);
        dPdu = (P -self.cache.P{t, k}) / epsilon;
        scale = 1;  % NOTE self.finite_diff() takes care of these
        return

      case 'eig'
        dPdu = dPdu_eig(self.cache.H_v{t, k}, self.cache.H_eig_factor{t, k}, self.system.B{k, c});

      case 'aux'
        dPdu = dPdu_auxmatrix(self.seq.tau(t) * self.cache.H{t, k}, self.system.B{k, c});

      case 'series_ss'
        %self.gradient_test(t, k, c);
        W = self.cache.W{t, k};
        n = length(W);  % number of squarings used for P
        s = 2^n;  % corresponding scale
        G = (self.seq.tau(t) / s) * self.cache.H{t, k};  % scaling

        % temp == expm(G)
        if n == 0
            temp = self.cache.P{t, k};
        else
            temp = W{1};
        end
        % d(expm(G))/du,  NOTE missing tau here
        dPdu = dPdu_series(G, self.system.B{k, c} / s, 7) * temp;
        % squaring
        for r=1:n
            dPdu = dPdu * W{r} +W{r} * dPdu;  % chain rule for the derivative
        end

      case 'series'
        % NOTE missing tau here
        dPdu = dPdu_series(self.seq.tau(t) * self.cache.H{t, k}, self.system.B{k, c}, 12);
        ut = t+1; % missing P{t} for dPdu
    end

    % NOTE every dPdu above except the fd one is missing a tau factor, which we add here
    scale = self.seq.tau(t) * self.seq.fields_deriv(t, c);
end
