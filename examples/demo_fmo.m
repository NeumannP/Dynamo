function dyn = demo_fmo(T)
% Example: Optimizing exciton transport in the FMO complex.

% Model from M.B. Plenio and S.F. Huelga, "Dephasing-assisted transport: quantum networks and biomolecules", NJP 10, 113019 (2008).
% Dipole moments etc. from F.Caruso et al., "Coherent open-loop optimal control of light-harvesting dynamics", arXiv:1103.0929v1.

% Ville Bergholm 2012-2014


SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;


%% Define the physics of the problem

% Exciton transport network with the fiducial sink as the last site,
% sucking probability from site 3

n_sites = 7;
dim = 2 * ones(1, n_sites + 1); % dimension vector: qubits, or spin-1/2 chain

% 0: ground, 1: excited state, SP lowers/annihilates
a = SP;
n_op = a' * a; % == (I-SZ) / 2; % number op


%% cleverness

% The Hamiltonian preserves exciton number, whereas the noise
% processes we use here leave the 0-or-1 exciton manifold
% invariant.

ddd  = n_sites + 2; % zero- and single-exciton subspaces (sink included)
mask = [1, 1+2.^(n_sites:-1:0)]; % states to keep: zero, single exciton at each site/sink
%away = setdiff(1:prod(dim), p); % states to throw away



%% parameters

qit.c = 299792458;       % speed of light in vacuum, m/s

% The unit for energy is 2 pi hbar c 100/m,
% and the unit for time (2 pi c 100/m)^{-1} \approx 5.31 ps.
TU = 1e12 / (2*pi * qit.c * 100) % ps

if nargin < 1
    % total sequence duration, in ps
    T = 5 * TU;
end

% dephasing
gamma = 0 * ones(1, n_sites)
% relaxation
Gamma = 0.5/188 * ones(1, n_sites)
% transfer rate from site 3 to sink
transfer_rate = 10/1.88


%% use a hybrid basis?

hybrid = 0
U = eye(ddd);
if hybrid
    basis = 'hybrid'
    if n_sites == 7
        ind = 1 + [1 2]; % rotate states 1 and 2
        R = R_y(-0.5*pi);
    else
        ind = 1 + [1 2 8];
        R = [-0.4261    0.6731    0.6045;
             -0.9001   -0.3827   -0.2083;
              0.0911   -0.6329    0.7689];
    end
    U(ind, ind) = R;
else
    basis = 'site';
end

desc = sprintf('FMO complex, %d sites, %s basis, transfer rate = %5.5g, Gamma = %5.5g',...
               n_sites, basis, transfer_rate, Gamma(1));


%% Lindblad ops

% NOTE the factors of sqrt(2) to match the parameters with our Lindblad eq. convention.

% dissipation and dephasing
diss = cell(1, n_sites);
deph = cell(1, n_sites);
for k = 1:n_sites
    diss{k} = op_list({{sqrt(2 * Gamma(k)) * a,    k}}, dim);
    diss{k} = clip(diss{k});

    deph{k} = op_list({{sqrt(2 * gamma(k)) * n_op, k}}, dim);
    deph{k} = clip(deph{k});
end

% sink
sink = {op_list({{a, 3; sqrt(2 * transfer_rate) * a', n_sites+1}}, dim)};
sink{1} = clip(sink{1});


%% drift Hamiltonian

if n_sites == 7
% Hamiltonian from Adolphs and Renger, Biophys. J. 91, 2778 (2006).
% Energies shifted by 12230 1/cm
H_FMO = [215, -104.1,  5.1,  -4.3,   4.7, -15.1,  -7.8;
           0,    220, 32.6,   7.1,   5.4,   8.3,   0.8;
           0,      0,    0, -46.8,   1.0,  -8.1,   5.1;
           0,      0,    0,   125, -70.7, -14.7, -61.5;
           0,      0,    0,     0,   450,  89.7,  -2.5;
           0,      0,    0,     0,     0,   330,  32.7;
           0,      0,    0,     0,     0,     0,   280];

st_labels = {'loss', '1', '2', '3*', '4', '5', '6', '7', 'sink'};

% pure state transfer, 1 to sink
initial = '10000000';
final   = '00000001';


% dipole moments of the chromophores in Debye
mu = ...
  [-3.081,  2.119, -1.669;...
   -3.481, -2.083, -0.190;...
   -0.819, -3.972, -0.331;...
   -3.390,  2.111, -1.080;...
   -3.196, -2.361,  0.792;...
   -0.621,  3.636,  1.882;...
   -1.619,  2.850, -2.584];

% times maximum electric field strength in (energy unit)/D
mu = mu * 15; %

elseif n_sites == 8
% 8-site Hamiltonian from Schmidt am Busch et al., J. Phys. Chem. Lett. 2, 93 (2011).
% Energies: \epsilon_border = 4 (in vivo setting)
% Energies shifted by 12176 1/cm
% NOTE The paper has a typo, missing a minus sign on H(1,2), which is fixed here.
H_FMO = [417, -94.8,   5.5,   -5.9,   7.1, -15.1, -12.2,  39.5;
           0,   254,  29.8,    7.6,   1.6,  13.1,   5.7,   7.9;
           0,     0,     0,  -58.9,  -1.2,  -9.3,   3.4,   1.4;
           0,     0,     0,    279, -64.1, -17.4, -62.3,  -1.6;
           0,     0,     0,      0,   420,  89.5,  -4.6,   4.4;
           0,     0,     0,      0,     0,   276,  35.1,  -9.1;
           0,     0,     0,      0,     0,     0,   293, -11.1;
           0,     0,     0,      0,     0,     0,     0,   471];

st_labels = {'loss', '1', '2', '3*', '4', '5', '6', '7', '8', 'sink'};

% pure state transfer, 8 to sink
initial = '000000010';
final   = '000000001';

else
    error('Wrong number of sites.')
end

H_FMO = H_FMO + triu(H_FMO, 1)'
H_FMO = blkdiag(0, H_FMO, 0); % add loss state and sink
H_drift = U' * H_FMO * U; % hybrid basis


%% drift Liouvillian (noise / dissipation)

L_drift = superop_lindblad(sink, H_drift) +superop_lindblad(diss) +superop_lindblad(deph);


%% controls

% Control Hamiltonians / Liouvillians
H_ctrl = {};
c_labels = {};

% transformed controls?
control_par = {};
par = [0, 500];

if 1
% dephasing controls
for k = 1:n_sites
    temp = op_list({{sqrt(2) * n_op, k}}, dim);
    temp = clip(temp);
    H_ctrl{end+1} = superop_lindblad({temp});
    c_labels{end+1} = sprintf('Dz_%d', k);
    control_par{end+1} = par;
end
control_type = char('m' + zeros(1, length(H_ctrl)));

else

% Laser field control ops (xr, xi, yr, yi, zr, zi)
for k = 1:3
    temp_r = 0;
    temp_i = 0;
    for s = 1:n_sites
        temp_r = temp_r -mu(s,k) * op_list({{SX, s}}, dim);
        temp_i = temp_i -mu(s,k) * op_list({{SY, s}}, dim);
    end
    H_ctrl{end+1} = temp_r(p,p);
    H_ctrl{end+1} = temp_i(p,p);
end
control_type = '......';
control_par = {};
c_labels = {'Xr', 'Xi', 'Yr', 'Yi', 'Zr', 'Zi'};
end


%% initial and final states

% pure state transfer
initial = to_state(initial, dim); initial = U' * initial(mask);
final   = to_state(final, dim);   final   = U' * final(mask);


%% set up Dynamo

dyn = dynamo('open state overlap', initial, final, L_drift, H_ctrl);
dyn.system.set_labels(desc, st_labels, c_labels);

% use the expensive-but-reliable gradient method
dyn.config.epsilon = 1e-4;
dyn.config.gradient_func = @gradient_g_finite_diff;


%% set up controls

% optimized dephasing rates for different values of T

% Plenio's
f7_5_plenio = [469.34, 5.36, 99.13, 5.55, 114.86, 1.88, 291.08];
f7_inf_plenio = [27.4, 26.84, 1.22, 87.12, 99.59, 232.76, 88.35];

% From Dynamo
f7_1_const = [ 4.4841  161.7343    0.0  112.7218  275.8279  500.0    0.0];
f7_2_const = [16.2512  159.4469    0.0  110.9683  214.1849  314.2133 6.6669];
f7_5_const = [22.4646  158.4658    0.0  110.3513  202.4489  161.6336 13.8490];
f7_1000_const = [6.8757  159.9530  0.0  113.4639  234.2341  500.0    0.0];

f8_1_const = [7.0533   166.5064   59.9556  169.3896  119.6819  189.6102  154.8109    0.0];
f8_2_const = [140.4333    0.0  230.2109    0.0    0.0    0.0    0.0  104.8807];
f8_5_const = [14.1791  162.1146   65.7091  170.4137   78.4960  179.4678  118.5341  0.0003];
f8_1000_const = [9.7400 164.6028  62.1406  170.8336   98.7679  186.5018  142.8440    0.0];


bins = 1
T = T / TU; % T originally in ps, now in TU
tau_par = T * [1, 0];
dyn.seq_init(bins, tau_par, control_type, control_par);

dyn.easy_control(1e-3 * f7_5_plenio, 0.0, 0, false);
%dyn.easy_control(1e-3 * rand(size(f7_5_const)), 0.0, 0, false);


%% now do the actual search

% "before" plot
figure(); dyn.plot_X(0.001);

dyn.ui_open();
dyn.search();

% "after" plot
figure(); dyn.plot_X(0.001);
return


function safe_x = clip(safe_x)
% Apply the subspace mask and a possible basis transformation U to an operator.
  safe_x  = U' * safe_x(mask,mask) * U;
end
end


function s = to_state(s, dim)
% Returns the state vector corresponding to the ket label string s.

  % calculate the linear index
  n = length(dim);
  s = s - '0';
  if any(s >= dim)
      error('Invalid basis ket.')
  end
  muls = fliplr(circshift(cumprod(fliplr(dim)), [0 1]));
  muls(end) = 1;
  ind = muls*s.';
  s = zeros(prod(dim), 1);
  s(ind+1) = 1; % MATLAB indexing starts at 1
end
