function dPdu = dPdu_eig(H_v, H_eig_factor, B)
% Compute the derivative dP_t/du_c using the eigendecomposition of H*tau.
% For efficiency tau has been factored into H when computing the
% eigenfactors, hence this function actually returns (dP/du_c) / tau.

    temp = H_v' * B * H_v;
    dPdu_eigenbasis = temp .* H_eig_factor; % note the elementwise multiplication .* here
    dPdu = H_v * dPdu_eigenbasis * H_v';    % to computational basis
end
