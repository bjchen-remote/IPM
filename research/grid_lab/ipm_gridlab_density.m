function [density,parts] = ipm_gridlab_density(coordinate,model)
%IPM_GRIDLAB_DENSITY Evaluate a positive-half analytic monitor density.
%   DENSITY = IPM_GRIDLAB_DENSITY(S,MODEL) evaluates MODEL for S >= 0.
%   The main 'asymmetric_power_front' family combines an asymmetric
%   algebraic core, a narrow left-front component, and a broad bridge.
%   'log_spline_teacher' is a nonparametric teacher used to identify where
%   a final analytic family is under-resolved.

validateattributes(coordinate,{'numeric'}, ...
    {'vector','real','finite','nonnegative'},mfilename,'coordinate');
if ~isstruct(model) || ~isscalar(model) || ~isfield(model,'kind')
    error('ipm:gridlab:DensityModel', ...
        'model must be a scalar structure with a kind field.');
end
wasRow = isrow(coordinate);
s = double(coordinate(:));
kind = validatestring(model.kind, ...
    {'asymmetric_power_front','log_spline_teacher'},mfilename,'model.kind');

switch kind
    case 'asymmetric_power_front'
        required_fields(model,{'floor','combinePower','contrastCap', ...
            'core','front','bridge'},'model');
        validateattributes(model.floor,{'numeric'}, ...
            {'scalar','real','finite','positive'},mfilename,'model.floor');
        validateattributes(model.combinePower,{'numeric'}, ...
            {'scalar','real','finite','>=',1},mfilename, ...
            'model.combinePower');
        validateattributes(model.contrastCap,{'numeric'}, ...
            {'scalar','real','positive'},mfilename,'model.contrastCap');

        core = asymmetric_component(s,model.core,'model.core');
        front = symmetric_component(s,model.front,'model.front');
        bridge = symmetric_component(s,model.bridge,'model.bridge');
        power = model.combinePower;
        focused = (core.^power+front.^power+bridge.^power).^(1/power);
        if isfinite(model.contrastCap)
            focused = model.contrastCap*tanh(focused/model.contrastCap);
        end
        densityColumn = model.floor+focused;
        parts = struct('core',restore(core,wasRow), ...
            'front',restore(front,wasRow), ...
            'bridge',restore(bridge,wasRow), ...
            'focused',restore(focused,wasRow));

    case 'log_spline_teacher'
        required_fields(model,{'knots','logDensity'},'model');
        knots = double(model.knots(:));
        logDensity = double(model.logDensity(:));
        validateattributes(knots,{'numeric'}, ...
            {'vector','real','finite','nonnegative','increasing'}, ...
            mfilename,'model.knots');
        validateattributes(logDensity,{'numeric'}, ...
            {'vector','real','finite','numel',numel(knots)}, ...
            mfilename,'model.logDensity');
        if numel(knots) < 4 || min(s) < knots(1) || max(s) > knots(end)
            error('ipm:gridlab:TeacherDomain', ...
                ['The teacher needs at least four knots and must cover the ' ...
                'requested coordinate interval.']);
        end
        densityColumn = exp(interp1(knots,logDensity,s,'pchip'));
        parts = struct('teacher',restore(densityColumn,wasRow));
end

if any(~isfinite(densityColumn)) || any(densityColumn <= 0)
    error('ipm:gridlab:DensityPositivity', ...
        'The monitor density must remain finite and strictly positive.');
end
density = restore(densityColumn,wasRow);
end

function values = asymmetric_component(s,component,name)
required_fields(component,{'center','strength','leftWidth','rightWidth', ...
    'leftShape','rightShape','leftTailPower','rightTailPower'},name);
validateattributes(component.center,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,[name '.center']);
validateattributes(component.strength,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,[name '.strength']);
for field = {'leftWidth','rightWidth','leftShape','rightShape', ...
        'leftTailPower','rightTailPower'}
    validateattributes(component.(field{1}),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,[name '.' field{1}]);
end
values = zeros(size(s));
left = s <= component.center;
z = (component.center-s(left))/component.leftWidth;
values(left) = component.strength*(1+z.^component.leftShape) ...
    .^(-component.leftTailPower/component.leftShape);
right = ~left;
z = (s(right)-component.center)/component.rightWidth;
values(right) = component.strength*(1+z.^component.rightShape) ...
    .^(-component.rightTailPower/component.rightShape);
end

function values = symmetric_component(s,component,name)
required_fields(component, ...
    {'center','strength','width','shape','tailPower'},name);
validateattributes(component.center,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,[name '.center']);
validateattributes(component.strength,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,[name '.strength']);
for field = {'width','shape','tailPower'}
    validateattributes(component.(field{1}),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,[name '.' field{1}]);
end
z = abs(s-component.center)/component.width;
values = component.strength*(1+z.^component.shape) ...
    .^(-component.tailPower/component.shape);
end

function required_fields(value,names,label)
if ~isstruct(value) || ~isscalar(value)
    error('ipm:gridlab:DensityComponent', ...
        '%s must be a scalar structure.',label);
end
missing = names(~isfield(value,names));
if ~isempty(missing)
    error('ipm:gridlab:DensityComponent', ...
        '%s is missing: %s.',label,strjoin(missing,', '));
end
end

function value = restore(value,wasRow)
if wasRow
    value = value.';
end
end
