function report = ipm_accellab_test_continuous_flat_events(outputRoot)
%IPM_ACCELLAB_TEST_CONTINUOUS_FLAT_EVENTS Explicit threshold-plateau guards.
% These adversarial C1 data use fixed linear jet maps, not PDE derivative
% matrices. They test invalid-event handling; maintained-FD accuracy is
% covered separately by test_continuous_inner and is not inferred here.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
x = -2:4; y = linspace(0,1,9)'; wall = [0,.9,1,.9,.9,.9,0];
jets = [1,.2,0,0,0,-.5,-1]; Dx = zeros(numel(x),numel(x)); Dx(:,3) = jets';
Dy = ipm.mesh.fdMatrix(y,1,7); omega = (1-y.^2)*wall;
wallCase = ipm_accellab_continuous_inner_rates(omega,zeros(size(omega)),x,y,Dx,Dy,[-2,4]);
assert(~wallCase.valid && any(strcmp(wallCase.invalidReasons, ...
    'threshold_flat_cell_intersects_connected_wall_component')));
x = -2:2; y = [0,.2,.4,.6,1]'; wall = [0,0,1,0,0]; trace = [1,.9,.9,.9,.5]';
Dx = zeros(numel(x),numel(x)); Dy = zeros(numel(y),numel(y)); Dy(end,1) = -1;
omega = trace*wall;
verticalCase = ipm_accellab_continuous_inner_rates(omega,zeros(size(omega)),x,y,Dx,Dy,[-2,2]);
assert(~verticalCase.valid && any(strcmp(verticalCase.invalidReasons, ...
    'threshold_flat_cell_intersects_connected_vertical_component')));
report = struct('status','passed','wallCase',wallCase,'verticalCase',verticalCase, ...
    'jetMaps','Fixed linear synthetic operators constructed to represent exact C1 plateau events; not maintained FD or a PDE experiment.', ...
    'rhsEvaluations',0,'pdeSteps',0,'poissonOperatorBuilds',0);
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['continuous_flat_events_',token]); mkdir(destination);
save(fullfile(destination,'report.mat'),'report');
fid = fopen(fullfile(destination,'report.json'),'w'); assert(fid>=0); cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
fprintf('CONTINUOUS_FLAT_EVENTS %s\n',destination);
end
