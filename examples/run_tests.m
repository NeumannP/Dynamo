% Test script for DYNAMO

d = demo_tasks('closed ket');
d = demo_tasks('closed ket phase');
d = demo_tasks('closed state');
%d = demo_tasks('closed state_partial');
d = demo_tasks('closed gate');
d = demo_tasks('closed gate phase');
d = demo_tasks('closed gate_partial');
d = demo_tasks('open state');
d = demo_tasks('open state overlap');
d = demo_tasks('open state_partial');
d = demo_tasks('open gate');
%d = demo_tasks('open gate_partial');
d = demo_tasks('abstract vector');
d = demo_tasks('abstract matrix');

disp('All tests passed.')
