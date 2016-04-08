function S = superop_comm(H)
% SUPEROP_COMM  Liouvillian superoperator equivalent for a commutator.
%  S = superop_comm(H);
%
%  [H, rho] == inv_vec(superop_comm(H) * vec(rho))

% Ville Bergholm 2011


S = lmul(H) -rmul(H);
