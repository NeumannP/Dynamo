function dyn = demo_coop_gates()
% Example: Concurrent optimization method of coop gates from
%! M. Braun and S. Glaser, New J. Phys. 16, 115002 (2014).

% Ville Bergholm 2015-2016


% Idea: Ramsey sequence: init to |0>, pi/2_y, wait T, -pi/2_y, measure |0><0|.
% The +-pi/2_y rotations do not have to be exact y rotations, it's
% enough if they bring you to the equator and back with extra
% phases that cancel each other. This extra (gauge?) freedom makes the gates
% shorter and/or easier to optimize.

global qit

% Liouvillian superoperator, projects to the equator of the Bloch sphere.
P_xy = [1, 0, 0, 1; 0, 2, 0, 0; 0, 0, 2, 0; 1, 0, 0, 1] / 2;

dim = 2;  % one qubit


%% Set up Dynamo

task = 'open state overlap';
ini = [1, 0].';
fin = ini;
H_drift = qit.sz / 2;
[H_ctrl, c_labels] = control_ops(dim, 'xy');


desc = task;
T = 2;
n_bins = 7;

dyn = dynamo(task, ini, fin, H_drift, H_ctrl);
dyn.system.set_labels(desc, dim, c_labels);
dyn.seq_init(n_bins, T * [0.5, 1.5]); %, control_type, control_par);

% random, constant initial controls
dyn.easy_control(0.1 * randn(1, 2));

mask = dyn.full_mask(false);

t = 4;  % P_xy slice
mask(t,:) = false;  % do not change the P_xy
dyn.cache.set_P(t, 1, P_xy);

dyn.ui_open();
dyn.search(mask);
