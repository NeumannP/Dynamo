function dPdu = dPdu_auxmatrix(Gtau, B)
% Compute the derivative dP/du_c using the auxiliary matrix method.
% For efficiency this function actually returns (dP/du_c) / tau.

    d = size(B);
    temp = expm([Gtau, B; zeros(d), Gtau]);
    dPdu = temp(1:d, d+1:end);
end
