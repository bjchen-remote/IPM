function [candidate,report] = ipm_accellab_secant(history,evaluate,user)
%IPM_ACCELLAB_SECANT Residual-tested reduced-space profile extrapolation.
%   Columns of HISTORY are same-coordinate states, oldest first. EVALUATE
%   returns [residual,audit], where audit.valid is a logical scalar. A fixed
%   linear residual metric must be used for every evaluation. This function
%   does not integrate time, create a solver state, or write a checkpoint.
%   User callbacks project and accept permit an exact gauge retraction and
%   independent residual gates. Every trial is freshly evaluated after its
%   projection; a failed trial is recorded and leaves the baseline intact.

if nargin < 3, user = struct(); end
opts = options(user);
validateattributes(history,{'double'}, ...
    {'2d','real','finite','nonempty'});
assert(size(history,2) >= 2,'ipm:AccelHistory', ...
    'At least two same-coordinate states are required.');
assert(isa(evaluate,'function_handle'),'ipm:AccelCallback', ...
    'evaluate must be a function handle.');
history = history(:,max(1,end-opts.memory):end);
count = size(history,2);
audits = cell(1,count);
[value,audits{1}] = evaluate(history(:,1));
validate_residual(value,audits{1});
residuals = zeros(numel(value),count);
residuals(:,1) = value(:);
for index = 2:count
    [value,audits{index}] = evaluate(history(:,index));
    validate_residual(value,audits{index});
    assert(numel(value) == size(residuals,1),'ipm:AccelResidualSize', ...
        'The residual dimension must remain fixed.');
    residuals(:,index) = value(:);
end
candidate = history(:,end);
baselineResidual = residuals(:,end);
baselineNorm = norm(baselineResidual);
report = struct('kind','stationary_profile_candidate_only', ...
    'accepted',false,'status','no_decreasing_trial', ...
    'isTimeTrajectory',false,'baselineNorm',baselineNorm, ...
    'candidateNorm',baselineNorm,'residualRatio',1, ...
    'rank',0,'singularValues',[],'coefficients',[], ...
    'coefficientNorm',0,'damping',0,'relativeJump',0, ...
    'evaluations',count,'baselineAudit',audits{end}, ...
    'candidateAudit',audits{end},'attempts',struct([]), ...
    'options',rmfield(opts,{'project','accept'}));
if baselineNorm == 0
    report.status = 'zero_baseline_residual';
    report.residualRatio = 0;
    return;
end
deltaState = diff(history,1,2);
deltaResidual = diff(residuals,1,2);
[left,sigma,right] = svd(deltaResidual,'econ');
singularValues = diag(sigma);
report.singularValues = singularValues;
if isempty(singularValues) || singularValues(1) == 0
    report.status = 'no_resolved_secant_direction';
    return;
end
keep = singularValues > opts.svdRelativeFloor*singularValues(1);
report.rank = nnz(keep);
if ~any(keep)
    report.status = 'no_resolved_secant_direction';
    return;
end
s = singularValues(keep);
filter = s./(s.^2+(opts.regularization*singularValues(1))^2);
coefficients = right(:,keep)*(filter.*(left(:,keep)'*baselineResidual));
report.coefficients = coefficients;
report.coefficientNorm = norm(coefficients);
if any(~isfinite(coefficients)) || norm(coefficients) > opts.maxCoefficientNorm
    report.status = 'coefficient_guard';
    return;
end
direction = -deltaState*coefficients;
for damping = opts.damping
    attempt = struct('damping',damping,'relativeJump',NaN, ...
        'residualNorm',NaN,'accepted',false,'reason','', ...
        'exceptionIdentifier','','exceptionMessage','','audit',struct());
    try
        trial = opts.project(history(:,end)+damping*direction);
        assert(isequal(size(trial),size(candidate)) && ...
            isreal(trial) && all(isfinite(trial)), ...
            'ipm:AccelProjection','Projection returned an invalid vector.');
        attempt.relativeJump = norm(trial-history(:,end))/ ...
            max(norm(history(:,end)),realmin);
        if attempt.relativeJump > opts.maxRelativeJump
            attempt.reason = 'jump_guard';
        else
            report.evaluations = report.evaluations+1;
            [trialResidual,trialAudit] = evaluate(trial);
            attempt.audit = trialAudit;
            validate_residual(trialResidual,trialAudit);
            assert(numel(trialResidual) == numel(baselineResidual), ...
                'ipm:AccelResidualSize','Trial changed residual dimension.');
            attempt.residualNorm = norm(trialResidual);
            sufficientDecrease = attempt.residualNorm <= ...
                (1-opts.minimumDecrease*damping)*baselineNorm;
            independentPassed = opts.accept(trialAudit,audits{end});
            assert(islogical(independentPassed) && isscalar(independentPassed), ...
                'ipm:AccelAccept','accept must return a logical scalar.');
            if sufficientDecrease && independentPassed
                candidate = trial;
                report.accepted = true;
                report.status = 'accepted_stationary_candidate';
                report.candidateNorm = attempt.residualNorm;
                report.residualRatio = attempt.residualNorm/baselineNorm;
                report.candidateAudit = trialAudit;
                report.damping = damping;
                report.relativeJump = attempt.relativeJump;
                attempt.accepted = true;
                attempt.reason = 'fresh_residual_decrease';
            else
                attempt.reason = 'residual_or_independent_guard';
            end
        end
    catch exception
        attempt.reason = 'trial_evaluation_rejected';
        attempt.exceptionIdentifier = exception.identifier;
        attempt.exceptionMessage = exception.message;
    end
    if isempty(report.attempts)
        report.attempts = attempt;
    else
        report.attempts(end+1) = attempt;
    end
    if report.accepted, return; end
end
end

function validate_residual(residual,audit)
assert(isnumeric(residual) && isreal(residual) && isvector(residual) && ...
    ~isempty(residual) && all(isfinite(residual(:))), ...
    'ipm:AccelResidual','Residual must be a finite real vector.');
assert(isstruct(audit) && isscalar(audit) && isfield(audit,'valid') && ...
    islogical(audit.valid) && isscalar(audit.valid) && audit.valid, ...
    'ipm:AccelAudit','Residual evaluation did not pass its validity audit.');
end

function opts = options(user)
opts = struct('memory',5,'svdRelativeFloor',1e-10, ...
    'regularization',1e-10,'maxCoefficientNorm',1e4, ...
    'maxRelativeJump',0.25,'minimumDecrease',1e-3, ...
    'damping',[1,0.5,0.25,0.125,0.0625], ...
    'project',@(x)x,'accept',@(trial,baseline)true);
assert(isstruct(user) && isscalar(user),'ipm:AccelOptions', ...
    'Options must be a scalar structure.');
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:AccelOptions', ...
    'An unknown acceleration option was supplied.');
for index = 1:numel(names), opts.(names{index}) = user.(names{index}); end
validateattributes(opts.memory,{'numeric'},{'scalar','integer','>=',1,'<=',20});
validateattributes(opts.svdRelativeFloor,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(opts.regularization,{'numeric'},{'scalar','finite','>=',0});
validateattributes(opts.maxCoefficientNorm,{'numeric'},{'scalar','finite','positive'});
validateattributes(opts.maxRelativeJump,{'numeric'},{'scalar','finite','positive'});
validateattributes(opts.minimumDecrease,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(opts.damping,{'numeric'},{'vector','nonempty','>',0,'<=',1});
opts.damping = opts.damping(:)';
assert(isa(opts.project,'function_handle') && isa(opts.accept,'function_handle'), ...
    'ipm:AccelOptions','project and accept must be function handles.');
end
