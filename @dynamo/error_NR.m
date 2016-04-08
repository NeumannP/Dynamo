function [L, J] = error_NR(self, control_mask)
% The Newton-Raphson (vector valued) error function L and its Jacobian.
% L is denoted by fraktur L in the docs.

% TODO instead of g, cache L? now there's no benefit in storing an
% intermediate value or is there?

% FIXME not functional
X_n = X_f' * self.X();
L = logm(X_n);
  
if nargin == 2
    J = self.gradient_NR(control_mask, L);
end

L = vec(project_to_su(L));
end
