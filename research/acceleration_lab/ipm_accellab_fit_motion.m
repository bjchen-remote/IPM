function [report,residual] = ipm_accellab_fit_motion(problem,user)
%IPM_ACCELLAB_FIT_MOTION Read-only weighted fit of inner coordinate motion.
%   Fits forcing + aPrime*translation + beta*scaleX + gamma*scaleY.
%   All columns, masks, weights, and amplitude normalization are fixed inputs.
if nargin < 2, user = struct(); end
opts = struct('fixedTranslationRate',NaN,'svdRelativeFloor',1e-10, ...
    'maximumCondition',1e8,'minimumFitPoints',12);
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:MotionOptions','Unknown motion-fit option.');
for index = 1:numel(names), opts.(names{index}) = user.(names{index}); end
validateattributes(opts.svdRelativeFloor,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(opts.maximumCondition,{'numeric'},{'scalar','finite','>=',1});
validateattributes(opts.minimumFitPoints,{'numeric'},{'scalar','integer','>=',3});
assert(isscalar(opts.fixedTranslationRate) && isreal(opts.fixedTranslationRate) && ...
    (isnan(opts.fixedTranslationRate) || isfinite(opts.fixedTranslationRate)), ...
    'ipm:MotionOptions','fixedTranslationRate must be finite or NaN.');
required = {'forcing','translation','scaleX','scaleY','weights','fitMask','masks','amplitudeScale'};
assert(isstruct(problem) && isscalar(problem) && all(isfield(problem,required)), ...
    'ipm:MotionProblem','Incomplete motion-fit problem.');
shape = size(problem.forcing);
for name = {'forcing','translation','scaleX','scaleY','weights'}
    value = problem.(name{1});
    assert(isnumeric(value) && isreal(value) && isequal(size(value),shape) && ...
        all(isfinite(value(:))),'ipm:MotionField','Invalid motion field.');
end
assert(all(problem.weights(:) > 0),'ipm:MotionWeights','Positive weights are required.');
validateattributes(problem.amplitudeScale,{'numeric'},{'scalar','finite','positive'});
mask = problem.fitMask;
assert(islogical(mask) && isequal(size(mask),shape) && ...
    nnz(mask) >= opts.minimumFitPoints,'ipm:MotionFitMask','Too few training nodes.');
columns = [problem.translation(mask),problem.scaleX(mask),problem.scaleY(mask)];
forcing = problem.forcing(mask);
fixed = ~isnan(opts.fixedTranslationRate);
if fixed
    forcing = forcing+opts.fixedTranslationRate*columns(:,1);
    columns = columns(:,2:3);
end
weight = sqrt(problem.weights(mask)/sum(problem.weights(mask)));
matrix = weight.*columns;
rightHandSide = -weight.*forcing;
columnNorms = vecnorm(matrix);
nonzero = all(columnNorms > realmin);
safeNorms = max(columnNorms,realmin);
[left,sigma,right] = svd(matrix./safeNorms,'econ');
singular = diag(sigma);
rankValue = nnz(singular > opts.svdRelativeFloor*max(singular));
if min(singular) <= realmin
    condition = Inf;
else
    condition = max(singular)/min(singular);
end
valid = nonzero && rankValue == size(columns,2) && condition <= opts.maximumCondition;
rates = NaN(3,1);
if valid
    reduced = (right*((left'*rightHandSide)./singular))./safeNorms';
    if fixed, rates = [opts.fixedTranslationRate;reduced]; else, rates = reduced; end
end
[metrics,residual] = ipm_accellab_motion_metrics(problem,rates);
report = struct('kind','read_only_inner_coordinate_motion_fit','valid',valid, ...
    'translationRate',rates(1),'logScaleXRate',rates(2),'logScaleYRate',rates(3), ...
    'fixedTranslation',fixed,'rank',rankValue,'scaledCondition',condition, ...
    'columnWeightedNorms',columnNorms,'scaledSingularValues',singular, ...
    'fitNodeCount',nnz(mask),'metrics',metrics,'options',opts, ...
    'isEvolutionGauge',false,'isTimeTrajectory',false);
end
