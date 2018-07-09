function J = angular_momentum(n)
% ANGULAR_MOMENTUM  Angular momentum matrices.
%  J = {Jx, Jy, Jz} = angular_momentum(d)
%
%  Returns the angular momentum matrices \vec(J)/\hbar
%  for the d-dimensional subspace defined by the
%  quantum number j == (d-1)/2, as a cell vector.
%
%  Normalization: [Jx, Jy] = 1i * Jz,  expm(-1i * J{k} * \theta) is a rotation by \theta radians.

% Ville Bergholm 2009-2012


if (n < 1)
  error('Dimension must be one or greater.')
end

j = (n-1)/2; % angular momentum quantum number, n == 2*j + 1

% raising operator in subspace J^2 = j*(j+1)
m = j;
Jplus = sparse(n,n);
for k=1:n-1
  m = m - 1;
  Jplus(k, k+1) = sqrt(j*(j+1) -m*(m+1));
end

% lowering operator
Jminus = Jplus';

% Jplus  = Jx + i*Jy
% Jminus = Jx - i*Jy

J = cell(1,3);
J{1} = 0.5*(Jplus + Jminus);
J{2} = 0.5*i*(Jminus - Jplus);
J{3} = spdiags([j:-1:-j].', 0, n, n);
