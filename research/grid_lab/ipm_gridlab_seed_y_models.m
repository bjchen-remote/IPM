function seeds = ipm_gridlab_seed_y_models(profiles,options)
%IPM_GRIDLAB_SEED_Y_MODELS Fit boundary-centred analytic y monitors.
%   SEEDS = IPM_GRIDLAB_SEED_Y_MODELS(PROFILES) uses the frozen positive-y
%   source profiles produced by IPM_GRIDLAB_FEATURE_PROFILES.  The analytic
%   candidates reuse the ASYMMETRIC_POWER_FRONT density family, but their
%   core is centred exactly at the y=0 boundary and the narrow front
%   component is disabled.  A broad bridge covers the union of the source
%   and envelope tails over all frozen times, so no interior left-front
%   geometry is assumed.
%
%   Both feature-profile schema versions 1 and 2 are accepted.  In either
%   case every record must carry yPositive, verticalSource,
%   verticalEnvelope, yPeakValue, and yCoreWidth.

if nargin < 2 || isempty(options)
    options = struct();
end
options = resolve_options(options);
records = validate_profiles(profiles);
recordCount = numel(records);

yCoordinates = cell(recordCount,1);
featureIndicators = cell(recordCount,1);
minimumFrozenDy = zeros(recordCount,1);
coreWidths = zeros(recordCount,1);
sourceSupportEdges = zeros(recordCount,1);
envelopeSupportEdges = zeros(recordCount,1);
peakValues = zeros(recordCount,1);
halfWidth = NaN;

for index = 1:recordCount
    [y,source,envelope,coreWidth,peakValue] = ...
        record_values(records(index),index);
    if index == 1
        halfWidth = y(end);
    elseif abs(y(end)-halfWidth) > ...
            1e3*eps(max(1,max(abs([y(end),halfWidth]))))
        error('ipm:gridlab:YSeedDomain', ...
            'All positive-y feature profiles must share one endpoint.');
    end

    sourceIndicator = normalize_signal(source);
    envelopeIndicator = normalize_signal(envelope);
    featureIndicator = combine_signals( ...
        sourceIndicator,envelopeIndicator,options);
    yCoordinates{index} = y;
    featureIndicators{index} = featureIndicator;
    minimumFrozenDy(index) = min(diff(y));
    coreWidths(index) = coreWidth;
    sourceSupportEdges(index) = support_edge( ...
        y,sourceIndicator,options.supportLevel);
    envelopeSupportEdges(index) = support_edge( ...
        y,envelopeIndicator,options.supportLevel);
    peakValues(index) = peakValue;
end

if ~any(cellfun(@(value) any(value > 0),featureIndicators))
    error('ipm:gridlab:YSeedSignal', ...
        'The positive-y source and envelope profiles are all zero.');
end
minimumWidth = max(options.minimumWidthFraction*halfWidth, ...
    options.minimumCellWidths*median(minimumFrozenDy));
if minimumWidth >= halfWidth
    error('ipm:gridlab:YSeedWidth', ...
        ['The requested minimum width leaves no positive-y interval for ' ...
        'fitting an analytic monitor.']);
end

% Use the widest observed half-maximum core.  This deliberately favours a
% broad boundary layer over a narrow monitor fitted to one frozen time.
coreWidth = max(max(coreWidths),minimumWidth);
minimumBridgeExtent = min(halfWidth,max(2*coreWidth,4*minimumWidth));
bridgeOuterEdge = min(halfWidth,max([sourceSupportEdges; ...
    envelopeSupportEdges;minimumBridgeExtent]));
bridgeCenter = 0.5*bridgeOuterEdge;
bridgeWidth = max(0.5*bridgeOuterEdge,minimumWidth);

base = struct();
base.name = 'y_boundary_power_bridge_seed';
base.kind = 'asymmetric_power_front';
base.floor = 1;
base.combinePower = options.combinePower;
base.contrastCap = options.contrastCaps(1);
base.core = struct('center',0,'strength',options.coreStrength, ...
    'leftWidth',coreWidth,'rightWidth',coreWidth, ...
    'leftShape',2,'rightShape',2, ...
    'leftTailPower',options.coreTailPower, ...
    'rightTailPower',options.coreTailPower);
% The density family requires a front component.  Its zero strength is an
% explicit statement that a boundary-centred y peak has no interior left
% front analogous to the positive-x wall-source profile.
base.front = struct('center',0,'strength',0, ...
    'width',minimumWidth,'shape',4,'tailPower',3);
base.bridge = struct('center',bridgeCenter, ...
    'strength',options.bridgeStrength,'width',bridgeWidth, ...
    'shape',2,'tailPower',options.bridgeTailPower);

candidateCount = numel(options.contrastCaps)* ...
    numel(options.coreWidthFactors)*numel(options.bridgeWidthFactors);
if candidateCount > 128
    error('ipm:gridlab:YSeedOptions', ...
        'The y-monitor option product exceeds 128 candidates.');
end
candidates = cell(candidateCount,1);
candidateIndex = 0;
for cap = options.contrastCaps
    for coreFactor = options.coreWidthFactors
        for bridgeFactor = options.bridgeWidthFactors
            candidateIndex = candidateIndex+1;
            model = base;
            model.contrastCap = cap;
            model.core.leftWidth = coreWidth*coreFactor;
            model.core.rightWidth = coreWidth*coreFactor;
            model.bridge.width = bridgeWidth*bridgeFactor;
            model.name = sprintf('y_boundary_cap_%g_core_%g_bridge_%g', ...
                cap,coreFactor,bridgeFactor);
            % Evaluate once here so malformed candidate parameters fail at
            % seed construction rather than inside a later grid sweep.
            ipm_gridlab_density(yCoordinates{1},model);
            candidates{candidateIndex} = model;
        end
    end
end

[teacher,reference] = build_teacher( ...
    yCoordinates,featureIndicators,options);
fitGeometry = struct( ...
    'boundaryCenter',0,'halfWidth',halfWidth, ...
    'minimumWidth',minimumWidth,'coreWidth',coreWidth, ...
    'coreWidthRange',[min(coreWidths),max(coreWidths)], ...
    'bridgeCenter',bridgeCenter,'bridgeWidth',bridgeWidth, ...
    'bridgeOuterEdge',bridgeOuterEdge, ...
    'sourceSupportRange', ...
        [min(sourceSupportEdges),max(sourceSupportEdges)], ...
    'envelopeSupportRange', ...
        [min(envelopeSupportEdges),max(envelopeSupportEdges)], ...
    'peakValueRange',[min(peakValues),max(peakValues)], ...
    'teacherPointCount',numel(reference), ...
    'frontComponentDisabled',true);
seeds = struct('schemaVersion',1, ...
    'kind','ipm_gridlab_y_seed_models', ...
    'profileSchemaVersion',profiles.schemaVersion, ...
    'coordinate','positive_y_boundary', ...
    'base',base,'teacher',teacher,'candidates',{candidates}, ...
    'fitGeometry',fitGeometry,'options',options);
end

function options = resolve_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:YSeedOptions', ...
        'options must be a scalar structure.');
end
defaults = struct( ...
    'minimumWidthFraction',1e-6,'minimumCellWidths',2, ...
    'combinePower',4,'maximumTeacherPoints',4097, ...
    'teacherContrast',12,'sourceWeight',1,'envelopeWeight',0.65, ...
    'supportLevel',0.1,'coreStrength',10,'bridgeStrength',4, ...
    'coreTailPower',2,'bridgeTailPower',2, ...
    'contrastCaps',[6,10,16], ...
    'coreWidthFactors',[1,1.5],'bridgeWidthFactors',1);
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(options),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:YSeedOptions', ...
        'Unknown option(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    end
end

validateattributes(options.minimumWidthFraction,{'numeric'}, ...
    {'scalar','real','finite','>',0,'<',1},mfilename, ...
    'options.minimumWidthFraction');
validateattributes(options.minimumCellWidths,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'options.minimumCellWidths');
validateattributes(options.combinePower,{'numeric'}, ...
    {'scalar','real','finite','>=',1},mfilename,'options.combinePower');
validateattributes(options.maximumTeacherPoints,{'numeric'}, ...
    {'scalar','integer','>=',17},mfilename, ...
    'options.maximumTeacherPoints');
validateattributes(options.supportLevel,{'numeric'}, ...
    {'scalar','real','finite','>',0,'<',1},mfilename, ...
    'options.supportLevel');
for field = {'teacherContrast','sourceWeight','envelopeWeight'}
    validateattributes(options.(field{1}),{'numeric'}, ...
        {'scalar','real','finite','nonnegative'},mfilename, ...
        ['options.' field{1}]);
end
if options.sourceWeight+options.envelopeWeight <= 0
    error('ipm:gridlab:YSeedOptions', ...
        'At least one vertical feature weight must be positive.');
end
for field = {'coreStrength','bridgeStrength','coreTailPower', ...
        'bridgeTailPower'}
    validateattributes(options.(field{1}),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename, ...
        ['options.' field{1}]);
end
for field = {'contrastCaps','coreWidthFactors','bridgeWidthFactors'}
    value = options.(field{1});
    validateattributes(value,{'numeric'}, ...
        {'vector','nonempty','real','finite','positive'},mfilename, ...
        ['options.' field{1}]);
    value = double(value(:)');
    if numel(unique(value)) ~= numel(value)
        error('ipm:gridlab:YSeedOptions', ...
            'options.%s must not contain duplicate values.',field{1});
    end
    options.(field{1}) = value;
end
if any(options.coreWidthFactors < 1) || ...
        any(options.bridgeWidthFactors < 1)
    error('ipm:gridlab:YSeedOptions', ...
        ['The y core and bridge width factors must be at least one; ' ...
        'this seed family does not narrow the observed boundary layer.']);
end
end

function records = validate_profiles(profiles)
if ~isstruct(profiles) || ~isscalar(profiles) || ...
        ~all(isfield(profiles,{'schemaVersion','kind','count','records'}))
    error('ipm:gridlab:YSeedProfiles', ...
        'Pass feature profiles with schema, kind, count, and records.');
end
validateattributes(profiles.schemaVersion,{'numeric'}, ...
    {'scalar','integer'},mfilename,'profiles.schemaVersion');
if ~ismember(profiles.schemaVersion,[1,2])
    error('ipm:gridlab:YSeedProfiles', ...
        'Feature-profile schemaVersion must be 1 or 2.');
end
if ~(ischar(profiles.kind) || ...
        (isstring(profiles.kind) && isscalar(profiles.kind))) || ...
        ~strcmp(char(profiles.kind),'ipm_frozen_grid_features')
    error('ipm:gridlab:YSeedProfiles', ...
        'profiles.kind must be ipm_frozen_grid_features.');
end
if ~isstruct(profiles.records) || isempty(profiles.records)
    error('ipm:gridlab:YSeedProfiles', ...
        'profiles.records must be a nonempty structure array.');
end
validateattributes(profiles.count,{'numeric'}, ...
    {'scalar','integer','positive'},mfilename,'profiles.count');
if profiles.count ~= numel(profiles.records)
    error('ipm:gridlab:YSeedProfiles', ...
        'profiles.count does not match the number of records.');
end
records = profiles.records(:);
required = {'yPositive','verticalSource','verticalEnvelope', ...
    'yPeakValue','yCoreWidth'};
missing = required(~isfield(records,required));
if ~isempty(missing)
    error('ipm:gridlab:YSeedProfiles', ...
        'Positive-y feature record(s) are missing: %s.', ...
        strjoin(missing,', '));
end
end

function [y,source,envelope,coreWidth,peakValue] = ...
        record_values(record,index)
validateattributes(record.yPositive,{'numeric'}, ...
    {'vector','real','finite','nonnegative','increasing'},mfilename, ...
    sprintf('profiles.records(%d).yPositive',index));
y = double(record.yPositive(:)');
if numel(y) < 4
    error('ipm:gridlab:YSeedRecord', ...
        'Each positive-y profile needs at least four points.');
end
tolerance = 1e3*eps(max(1,max(abs(y))));
if abs(y(1)) > tolerance || y(end) <= 0
    error('ipm:gridlab:YSeedDomain', ...
        'Each yPositive vector must start at zero and end above zero.');
end
y(1) = 0;
validateattributes(record.verticalSource,{'numeric'}, ...
    {'vector','real','finite','nonnegative','numel',numel(y)}, ...
    mfilename,sprintf('profiles.records(%d).verticalSource',index));
validateattributes(record.verticalEnvelope,{'numeric'}, ...
    {'vector','real','finite','nonnegative','numel',numel(y)}, ...
    mfilename,sprintf('profiles.records(%d).verticalEnvelope',index));
source = double(record.verticalSource(:)');
envelope = double(record.verticalEnvelope(:)');
validateattributes(record.yPeakValue,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename, ...
    sprintf('profiles.records(%d).yPeakValue',index));
validateattributes(record.yCoreWidth,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    sprintf('profiles.records(%d).yCoreWidth',index));
if record.yCoreWidth > y(end)+tolerance
    error('ipm:gridlab:YSeedRecord', ...
        'A yCoreWidth cannot exceed its positive-y domain.');
end
coreWidth = double(record.yCoreWidth);
peakValue = double(record.yPeakValue);
end

function values = normalize_signal(values)
scale = max(values);
if scale <= 100*eps(max(1,norm(values,inf)))
    values = zeros(size(values));
else
    values = values/scale;
end
end

function combined = combine_signals(source,envelope,options)
power = options.combinePower;
combined = ((options.sourceWeight*source).^power + ...
    (options.envelopeWeight*envelope).^power).^(1/power);
combined = normalize_signal(combined);
end

function edge = support_edge(y,signal,level)
active = find(signal >= level,1,'last');
if isempty(active)
    edge = 0;
elseif active == numel(y)
    edge = y(end);
else
    y1 = y(active);
    y2 = y(active+1);
    value1 = signal(active);
    value2 = signal(active+1);
    if value2 == value1
        edge = 0.5*(y1+y2);
    else
        fraction = (level-value1)/(value2-value1);
        fraction = min(max(fraction,0),1);
        edge = y1+fraction*(y2-y1);
    end
end
end

function [teacher,reference] = build_teacher( ...
        coordinates,indicators,options)
reference = unique([coordinates{:}],'sorted');
if numel(reference) > options.maximumTeacherPoints
    retained = unique(round(linspace( ...
        1,numel(reference),options.maximumTeacherPoints)));
    reference = reference(retained);
end
teacherStack = zeros(numel(coordinates),numel(reference));
for index = 1:numel(coordinates)
    teacherStack(index,:) = max(0,interp1( ...
        coordinates{index},indicators{index},reference,'pchip'));
end
teacherFeature = normalize_signal(max(teacherStack,[],1));
teacherDensity = 1+options.teacherContrast*teacherFeature;
teacher = struct('name','worst_time_y_log_spline_teacher', ...
    'kind','log_spline_teacher','knots',reference, ...
    'logDensity',log(teacherDensity));
end
