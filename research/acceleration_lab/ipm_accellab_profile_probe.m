function [rhoCandidate,report] = ipm_accellab_profile_probe(states,user)
%IPM_ACCELLAB_PROFILE_PROBE Try a stationary profile using the unchanged RHS.
%   STATES is a cell of same-grid schema-4 isotropic accepted solver states,
%   oldest first. Only quadratic-peak / transport-anchor is registered here.
%   The returned rho is a research candidate, never a continuation state.
%   No clock, state, config, ops, result, or checkpoint is mutated or saved.

if nargin < 2, user = struct(); end
opts = struct('coreX',[0,4],'coreY',[0,2], ...
    'gaugeTolerance',1e-8,'historyGaugeTolerance',1e-5, ...
    'secant',struct());
assert(isstruct(user) && isscalar(user),'ipm:AccelProfileOptions', ...
    'Options must be a scalar structure.');
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:AccelProfileOptions', ...
    'Unknown profile-probe option.');
for index = 1:numel(names), opts.(names{index}) = user.(names{index}); end
validateattributes(opts.coreX,{'numeric'},{'vector','numel',2,'finite','increasing'});
validateattributes(opts.coreY,{'numeric'},{'vector','numel',2,'finite','increasing'});
validateattributes(opts.gaugeTolerance,{'numeric'},{'scalar','finite','positive'});
validateattributes(opts.historyGaugeTolerance,{'numeric'},{'scalar','finite','positive'});
assert(iscell(states) && numel(states) >= 2 && numel(states) <= 21, ...
    'ipm:AccelProfileHistory','Supply between 2 and 21 solver states.');
base = states{end};
validate_state(base);
ops = base.ops;
for index = 1:numel(states)
    item = states{index};
    validate_state(item);
    groups = {'grid','physics','elliptic','transport','scaling'};
    for group = groups
        assert(isequaln(item.config.(group{1}),base.config.(group{1})), ...
            'ipm:AccelProfileContract','History changed frozen numerical data.');
    end
    assert(isequal(item.ops.x,ops.x) && isequal(item.ops.y,ops.y) && ...
        isequaln(item.ops.rescaling,ops.rescaling), ...
        'ipm:AccelProfileGrid','History changed its grid or exact gauge references.');
    if index > 1
        assert(item.scale.canonicalTime > states{index-1}.scale.canonicalTime, ...
            'ipm:AccelProfileTime','History canonical times must increase.');
    end
end
shape = size(base.rho);
targetPeak = peak_value(base.rho,ops);
weights = ops.integrationWeights;
assert(all(isfinite(weights(:))) && all(weights(:) > 0), ...
    'ipm:AccelProfileWeights','Positive quadrature weights are required.');
core = ops.X >= opts.coreX(1) & ops.X <= opts.coreX(2) & ...
    ops.Y >= opts.coreY(1) & ops.Y <= opts.coreY(2);
wall = ops.x >= opts.coreX(1) & ops.x <= opts.coreX(2);
assert(nnz(core) >= 9 && nnz(wall) >= 3, ...
    'ipm:AccelProfileCore','The fixed core window needs sufficient nodes.');
weights = weights/sum(weights,'all');
coreWeights = weights(core)/sum(weights(core));
wallWeights = ops.integrationWeights(1,wall);
wallWeights = wallWeights/sum(wallWeights);
baseX = base.rho*ops.Dx';
baseY = ops.Dy*base.rho;
densityScale = max(sqrt(sum(weights.*base.rho.^2,'all')),realmin);
gradientScale = max(sqrt(sum(coreWeights.* ...
    (baseX(core).^2+baseY(core).^2))),realmin);
wallScale = max(sqrt(sum(wallWeights.*baseX(1,wall).^2)),realmin);
history = zeros(numel(base.rho),numel(states));
for index = 1:numel(states)
    value = states{index}.rho;
    assert(abs(peak_value(value,ops)/targetPeak-1) <= ...
        opts.historyGaugeTolerance,'ipm:AccelProfileGaugeDrift', ...
        'History does not lie on the same resolved amplitude gauge surface.');
    history(:,index) = value(:);
end
secantOptions = opts.secant;
assert(~isfield(secantOptions,'project') && ~isfield(secantOptions,'accept'), ...
    'ipm:AccelProfileOptions','Profile projection and audit cannot be overridden.');
secantOptions.project = @project;
secantOptions.accept = @independent_gates;
[candidate,report] = ipm_accellab_secant(history,@evaluate,secantOptions);
rhoCandidate = reshape(candidate,shape);
report.profileOptions = opts;
report.scalingContract = base.config.scaling.scalingContract;
report.targetQuadraticPeak = targetPeak;
report.sourceCanonicalTimes = cellfun(@(s)s.scale.canonicalTime,states);
report.shapeDiagnostics = ipm_accellab_shape_history(states);
report.artifactDisposition = 'research_profile_candidate_no_time_or_checkpoint';
report.requiredNextEvidence = { ...
    'same discrete RHS residual reduction on independent holdout states', ...
    'short unmodified dynamics relaxation of a separately registered branch', ...
    'same physical profile under matched grid and box refinement'};

    function output = project(input)
        field = reshape(input,shape);
        amplitude = targetPeak/peak_value(field,ops);
        assert(isfinite(amplitude) && amplitude > 0, ...
            'ipm:AccelProfileProjection','Positive exact amplitude retraction failed.');
        output = amplitude*input;
    end

    function [vector,audit] = evaluate(input)
        field = reshape(input,shape);
        [rhs,flow] = ipm.evolve.flow(field,ops,base.scale);
        rhsX = rhs*ops.Dx';
        rhsY = ops.Dy*rhs;
        fullPart = sqrt(weights).*rhs/densityScale;
        coreXPart = sqrt(coreWeights).*rhsX(core)/gradientScale;
        coreYPart = sqrt(coreWeights).*rhsY(core)/gradientScale;
        wallPart = sqrt(wallWeights).*rhsX(1,wall)/wallScale;
        vector = [fullPart(:);coreXPart(:);coreYPart(:);wallPart(:)];
        audit = struct('valid',true, ...
            'fullRhsRelativeL2',norm(fullPart(:)), ...
            'coreGradientRhsRelativeL2',hypot(norm(coreXPart),norm(coreYPart)), ...
            'wallGradientRhsRelativeL2',norm(wallPart), ...
            'fullRhsRelativeInf',max(abs(rhs),[],'all')/densityScale, ...
            'coreGradientRhsRelativeInf',max(hypot(rhsX(core),rhsY(core)))/gradientScale, ...
            'cl',flow.canonicalCL,'cOmega',flow.canonicalCOmega, ...
            'quadraticPeak',flow.omegaGaugeQuadraticPeakValue, ...
            'omegaGaugeResidual',flow.omegaGaugeResidual, ...
            'lengthGaugeResidual',flow.lengthGaugeResidual, ...
            'poissonResidual',flow.poissonResidual);
        audit.valid = all(isfinite(vector)) && isfinite(flow.poissonResidual) && ...
            abs(flow.omegaGaugeResidual) <= opts.gaugeTolerance && ...
            abs(flow.lengthGaugeResidual) <= opts.gaugeTolerance && ...
            flow.timeSpeed == 1 && flow.widthCorrection == 0;
    end

    function passed = independent_gates(trial,baseline)
        fields = {'fullRhsRelativeL2','coreGradientRhsRelativeL2', ...
            'wallGradientRhsRelativeL2','fullRhsRelativeInf', ...
            'coreGradientRhsRelativeInf'};
        passed = trial.valid && ...
            abs(trial.quadraticPeak/targetPeak-1) <= opts.gaugeTolerance;
        for field = fields
            passed = passed && trial.(field{1}) <= ...
                baseline.(field{1})*(1+1e-12)+1e-14;
        end
    end
end

function validate_state(state)
required = {'rho','ops','config','scale'};
assert(isstruct(state) && isscalar(state) && all(isfield(state,required)), ...
    'ipm:AccelProfileState','A complete solver state is required.');
config = state.config;
s = config.scaling;
assert(config.schemaVersion == 4 && config.frozen && ...
    strcmp(s.scalingContract,'exact_gauge_no_feedback_v1') && ...
    strcmp(s.rescalingMode,'dynamic') && ...
    strcmp(s.dynamicScaleGeometry,'isotropic') && ...
    strcmp(s.lengthGauge,'transport_anchor') && s.transportAnchorX == 1 && ...
    strcmp(s.cOmegaGauge,'wall_omega_quadratic_peak') && ...
    s.lengthScaleGain == 1 && s.widthGaugeGain == 0 && ...
    s.omegaGaugeGain == 0 && s.travelingWaveGain == 0 && ...
    ~config.remesh.adaptiveRemesh && state.ops.remeshCount == 0, ...
    'ipm:AccelProfileContract','The registered exact fixed-grid gauge is required.');
assert(isequal(size(state.rho),[numel(state.ops.y),numel(state.ops.x)]) && ...
    isreal(state.rho) && all(isfinite(state.rho(:))), ...
    'ipm:AccelProfileState','State field and paired grid are inconsistent.');
end

function value = peak_value(rho,ops)
source = rho*ops.Dx';
tracking = abs(ops.x-ops.rescaling.pinX) <= ops.rescaling.peakTrackingHalfWidth;
indices = find(tracking);
assert(~isempty(indices),'ipm:AccelProfilePeak','Empty tracking window.');
[~,localIndex] = max(source(1,tracking));
peak = ipm.evolve.quadraticPeakFunctional(source(1,:),ops.x,indices(localIndex));
value = peak.value;
end
