function [H] = op_list(G, dim)
% OP_LIST  Operator consisting of k-local terms, given as a list.
%  O = op_list(G, dim)
%
%  Returns the operator O defined by the connection list G.
%  dim is a vector of subsystem dimensions for O.
%
%  G is a cell vector of cell arrays, G = {c_1, c_2, ..., c_n},
%  where each cell array c_i corresponds to a term in O.
%
%  An array that has 2 columns and k rows, c_i = {A1, s1; A2, s2; ... ; Ak, sk},
%  where Aj are operators and sj subsystem indices, corresponds to the
%  k-local term given by the tensor product
%
%    A1_{s1} * A2_{s2} * ... * Ak_{sk}.
%
%  The dimensions of all operators acting on subsystem sj must match dim(sj).
%
%  Alternatively one can think of G as defining a hypergraph, where
%  each subsystem corresponds to a vertex and each array c_i in the list
%  describes a hyperedge connecting the vertices {s1, s2, ..., sk}.
%
%  Example: The connection list
%    G = {{sz,1}, {sx,1;sx,3}, {sy,1;sy,3}, {sz,1;sz,3}, ...
%         {sz,2}, {A,2;B+C,3}, {2*sz,3}}
%
%  corresponds to the operator
%    O = sz_1 +sz_2 +2*sz_3 +sx_1*sx_3 +sy_1*sy_3 +sz_1*sz_3 +A_2*(B+C)_3.

% Ville Bergholm 2009-2017


if nargin < 2
  % TODO we could try to infer dim from the operators
  error('Need both the list and the dimension vector.')
end

n = length(G(:)); % number of terms
D = prod(dim);
H = sparse(D, D);

for k=1:n
  spec = G{k}; % m*2 cell array
  if isempty(spec)
      continue  % for convenience an empty spec is skipped
  end
  s = size(spec);
  if (s(2) ~= 2 || s(1) < 1)
    error('Malformed term spec %d', k)
  end

  a = 0; % last subsystem taken care of
  term = 1;
  for j=1:s(1)
    b = spec{j,2}; % subsystem number
    if (b <= a)
      error('Spec %d not in ascending order.', k)
    end
    if (size(spec{j,1}, 2) ~= dim(b))
      error('The dimension of operator %d in spec %d does not match dim.',j,k)
    end
    term = mkron(term, speye(prod(dim(a+1:b-1))), spec{j,1});
    a = b;
  end
  term = mkron(term, speye(prod(dim(a+1:end))));

  H = H + term;
end
