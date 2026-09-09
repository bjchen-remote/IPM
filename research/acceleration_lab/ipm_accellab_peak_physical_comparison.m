function report = ipm_accellab_peak_physical_comparison(integrationDirectory)
%IPM_ACCELLAB_PEAK_PHYSICAL_COMPARISON Common physical time, tiny fields only.
% Uses two fresh maintained spatial RHS values per saved trajectory, with
% its actual equation identity, for cubic temporal Hermite interpolation.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
[~,token] = fileparts(tempname); destination = fullfile(integrationDirectory,['physical_comparison_',token]);
mkdir(destination);
registration = struct('sourceDirectory',integrationDirectory,'outputDirectory',destination, ...
    'physicalTime',0.075,'physicalX',linspace(-3,3,97),'physicalY',linspace(0,3,49), ...
    'temporalInterpolator','cubic Hermite in physical time using genuine endpoint rho_tau / t_tau and scale rates / t_tau', ...
    'spatialInterpolator','linear-in-data C1 tensor cubic Hermite, with maintained Dx and Dy', ...
    'secondaryTemporalObserver','piecewise linear in physical time, reported without selecting between them', ...
    'physicalDensityFormula','R(Cl*x+Xshift,Cl*y)/Comega', ...
    'physicalGradientFormula','Cl/Comega times the maintained-Dx wall/bulk gradient interpolant', ...
    'pdeSteps',0,'scope','Common-time physical observations of existing tiny trajectories; two grids are not a continuum convergence certificate.');
save(fullfile(destination,'registration.mat'),'registration'); write_json(fullfile(destination,'registration.json'),registration);
[PX,PY] = meshgrid(registration.physicalX,registration.physicalY);
weights = trapezoid_weights(registration.physicalY)'*trapezoid_weights(registration.physicalX);
modes = {'native_quadratic','research_hermite_C1'}; sizes = [49,65];
observations = cell(2,2,3); entries = cell(2,2,3);
try
    for g = 1:2
        first = load(fullfile(integrationDirectory,sprintf('native_quadratic_n%d_level1.mat',sizes(g))),'run');
        opts = flat_options(first.run.operatorSourceConfig); state = ipm.evolve.initialize(opts);
        for m = 1:2
            for level = 1:3
                file = fullfile(integrationDirectory,sprintf('%s_n%d_level%d.mat',modes{m},sizes(g),level));
                data = load(file,'run'); run = data.run;
                assert(isequal(run.x,state.ops.x) && isequal(run.y,state.ops.y) && ...
                    isequal(run.Dx,state.ops.Dx) && isequal(run.Dy,state.ops.Dy));
                a = [run.accepted{:}]; times = [a.physicalTime]; target = registration.physicalTime;
                j = find(times <= target,1,'last');
                assert(~isempty(j) && j < numel(times) && all(diff(times) > 0));
                left = a(j); right = a(j+1); zl = pack(left); zr = pack(right);
                ql = run.rhoSnapshots{j}; qr = run.rhoSnapshots{j+1};
                el = ipm_accellab_peak_rhs(ql,zl,state.ops,run.mode);
                er = ipm_accellab_peak_rhs(qr,zr,state.ops,run.mode);
                assert(abs(el.diagnostic.H-left.H) < 1e-13 && abs(er.diagnostic.H-right.H) < 1e-13);
                width = times(j+1)-times(j); fraction = (target-times(j))/width;
                hl = exp(zl(2)-zl(1)); hr = exp(zr(2)-zr(1));
                rho = hermite(ql,qr,el.rhoRate/hl,er.rhoRate/hr,width,fraction);
                z = hermite(zl,zr,el.scaleRate/hl,er.scaleRate/hr,width,fraction);
                % This is observer interpolation at a known physical time,
                % not a state reset or a step in either trajectory.
                assert(abs(z(3)-target) < 1e-14);
                linearRho = (1-fraction)*ql+fraction*qr;
                linearZ = (1-fraction)*zl+fraction*zr;
                value = physical_observer(rho,z,state.ops,PX,PY);
                linear = physical_observer(linearRho,linearZ,state.ops,PX,PY);
                observations{g,m,level} = value;
                entries{g,m,level} = struct('mode',run.mode,'nx',sizes(g),'dt',run.dt, ...
                    'sourceFile',file,'physicalTime',target,'bracketRows',[j,j+1], ...
                    'physicalBracket',times(j:j+1),'fraction',fraction,'interpolatedScale',z, ...
                    'temporalLinearDifference',difference(linear,value,weights), ...
                    'densityDerivativeVsMaintainedGradientDifference',norms(value.densityDerivativeX-value.gradient,weights), ...
                    'newRhsEvaluations',2,'pdeSteps',0);
            end
        end
        clear state;
    end
    gaugeDifferences = cell(2,3); spatialDifferences = cell(2,3); timeDifferences = cell(2,2,2);
    for g = 1:2
        for level = 1:3
            gaugeDifferences{g,level} = difference(observations{g,2,level},observations{g,1,level},weights);
        end
    end
    for m = 1:2
        for level = 1:3
            spatialDifferences{m,level} = difference(observations{1,m,level},observations{2,m,level},weights);
        end
    end
    for g = 1:2
        for m = 1:2
            for level = 1:2
                timeDifferences{g,m,level} = difference(observations{g,m,level},observations{g,m,level+1},weights);
            end
        end
    end
    report = struct('status','completed_common_physical_time_observer','registration',registration, ...
        'entries',{entries},'gaugeDifferencesC1MinusNative',{gaugeDifferences}, ...
        'spatialDifferences49Minus65',{spatialDifferences},'adjacentTimeLevelDifferences',{timeDifferences}, ...
        'cellDimensionOrder',struct('entries','grid, mode, dtLevel','gaugeDifferences','grid, dtLevel', ...
        'spatialDifferences','mode, dtLevel','timeDifferences','grid, mode, adjacentDtPair'), ...
        'gridOrder',sizes,'modeOrder',{modes},'newBracketRhsEvaluations',24,'tinyOperatorInitializations',2, ...
        'newPdeSteps',0,'interpretation','Gauge-coordinate field error and physical field error are different. These fixed observations retain temporal interpolation sensitivity and finite-grid spatial differences; none certifies the continuum PDE.');
    save(fullfile(destination,'physical_fields.mat'),'observations','PX','PY','weights','-v7.3');
    save(fullfile(destination,'report.mat'),'report'); write_json(fullfile(destination,'report.json'),report);
catch exception
    failure = struct('registration',registration,'identifier',exception.identifier,'message',exception.message);
    save(fullfile(destination,'failure.mat'),'failure'); write_json(fullfile(destination,'failure.json'),failure); rethrow(exception);
end
fprintf('HERMITE_PEAK_PHYSICAL_COMPARISON %s\n',destination);
end

function value = physical_observer(rho,z,ops,PX,PY)
cl = exp(z(1)); cw = exp(z(2)); X = cl*PX+z(4); Y = cl*PY;
assert(min(X,[],'all') >= ops.x(1) && max(X,[],'all') <= ops.x(end) && ...
    min(Y,[],'all') >= ops.y(1) && max(Y,[],'all') <= ops.y(end));
a = ipm_accellab_tensor_hermite(rho,ops.x,ops.y,ops.Dx,ops.Dy,X,Y);
b = ipm_accellab_tensor_hermite(rho*ops.Dx',ops.x,ops.y,ops.Dx,ops.Dy,X,Y);
value = struct('density',a.value/cw,'gradient',cl/cw*b.value, ...
    'densityDerivativeX',cl/cw*a.derivativeX,'z',z);
end

function v = hermite(a,b,da,db,h,s)
v = (2*s^3-3*s^2+1)*a+(s^3-2*s^2+s)*h*da+ ...
    (-2*s^3+3*s^2)*b+(s^3-s^2)*h*db;
end

function d = difference(a,b,w)
d = struct('densityAbsoluteL2Infinity',norms(a.density-b.density,w), ...
    'densityRelativeL2Infinity',norms(a.density-b.density,w)./norms(b.density,w), ...
    'gradientAbsoluteL2Infinity',norms(a.gradient-b.gradient,w), ...
    'gradientRelativeL2Infinity',norms(a.gradient-b.gradient,w)./norms(b.gradient,w), ...
    'densityDerivativeAbsoluteL2Infinity',norms(a.densityDerivativeX-b.densityDerivativeX,w));
end

function n = norms(a,w)
n = [sqrt(sum(a.^2.*w,'all')/sum(w,'all')),max(abs(a),[],'all')];
end

function w = trapezoid_weights(x)
x = x(:)'; h = diff(x); w = [h(1)/2,(h(1:end-1)+h(2:end))/2,h(end)/2];
end

function z = pack(d)
z = [d.logC_l;d.logC_omega;d.physicalTime;d.X_shift;d.canonicalTime];
end

function opts = flat_options(config)
schema = ipm.config.schema(); opts = struct();
for k = 1:numel(schema.domainNames)
    domain = config.(schema.domainNames{k}); names = fieldnames(domain);
    for j = 1:numel(names), opts.(names{j}) = domain.(names{j}); end
end
end

function write_json(file,value)
fid = fopen(file,'w'); assert(fid >= 0); cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
