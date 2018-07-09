function seq = import(self, A, in_polar, TU)
% Imports a control sequence array into a Dynamo instance.
%    
%  The rows of A correspond to pulses, each row consists of
%    [tau (s), amp_1 (Hz), amp_2 (Hz), ...]
%  or, if in_polar is true,
%    [tau (s), f_rabi_1 (Hz), phi_1, f_rabi_2 (Hz), phi_2, ...]
%
%  If TU is given, the units for the data in A are TU s and 1/(TU s) instead.
%
%  Importing an exported sequence should leave Dynamo in the exact same state.

% Ville Bergholm 2013-2016


if nargin < 4
    % default: no explicit time unit for the input data
    TU = 1;
    if nargin < 3
        % default: import the controls as they are
        in_polar = false;
    end
end

% our internal time unit (in s)
TU_self = self.system.TU;
if isempty(TU_self)
    % default
    TU_self = 1;
end
scale = TU / TU_self;

% split the array up
tau = A(:, 1) * scale;
T = sum(tau)
A = A(:, 2:end);
n_bins = size(A, 1)
n_controls = size(A, 2)

if self.system.n_controls() ~= n_controls
    error('Wrong number of controls.')
end

if in_polar
    % transform the controls from polar to cartesian coordinates
    f = zeros(size(A));
    for k=1:2:n_controls
        amp = A(:, k) * 2*pi / scale;
        phi = A(:, k+1);
        f(:, k)   = amp .* cos(phi);
        f(:, k+1) = amp .* sin(phi);
    end
else
    % assume cartesian
    f = A * 2*pi / scale;
end

%self.seq_init(n_bins, T * [0.5, 1]);
self.set_controls(f, tau);
end
