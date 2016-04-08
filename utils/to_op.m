function psi = to_op(psi)
% TO_OP  Converts state vectors to state operators.
%  rho = to_op(psi)

% Ville Bergholm 2014

if size(psi, 2) == 1
    psi = psi * psi';
end
