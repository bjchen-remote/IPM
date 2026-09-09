function opts = strong_mesh_case(name)
%STRONG_MESH_CASE Reproducible options for accepted physical-mesh runs.
%   Pass the returned structure to the only simulation entry point:
%       result = ipm.solve(strong_mesh_case('production'));

if nargin < 1
    name = 'production';
end
common = struct('xlim',[-400,400],'ymax',400, ...
    'physicalFinalTime',1.44,'remeshMaximumCellRatio',6, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'saveResults',true,'verbose',false);
switch lower(char(name))
    case 'production'
        opts = common;
        opts.nx = 1025;
        opts.ny = 513;
        opts.targetCenterSpacing = [0.01,0.01];
        opts.maxDt = 2.5e-4;
        opts.resultFile = fullfile('result','production', ...
            'ipm_physical_strong_mesh_1025.mat');
    case 'spatial_513'
        opts = common;
        opts.nx = 513;
        opts.ny = 257;
        opts.targetCenterSpacing = [0.02,0.02];
        opts.maxDt = 2.5e-4;
        opts.physicalFinalTime = 1.43;
        opts.resultFile = fullfile('result','verification', ...
            'mesh_v2_strong_dt25_513.mat');
    otherwise
        error('ipm:UnknownStrongMeshCase', ...
            'Use ''production'' or ''spatial_513''.');
end
end
