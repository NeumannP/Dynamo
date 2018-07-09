function dPdu = dPdu_series(G, B, n)
% Compute the derivative dP/du using the commutator series method.
% Here P = expm(G), G = A +u*B.
% For efficiency this function actually returns (dP/du) * P^{-1}.

    dPdu = B;
    for k=2:n
        % commutator
        B = (G*B -B*G) / k;
        dPdu = dPdu +B;
    end
end
