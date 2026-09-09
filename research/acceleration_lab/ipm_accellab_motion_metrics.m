function [metrics,residual] = ipm_accellab_motion_metrics(problem,rates)
%IPM_ACCELLAB_MOTION_METRICS Apply prescribed diagnostic coordinate rates.
validateattributes(rates,{'numeric'},{'vector','numel',3,'real'});
residual = problem.forcing+rates(1)*problem.translation+ ...
    rates(2)*problem.scaleX+rates(3)*problem.scaleY;
metrics = struct();
maskNames = fieldnames(problem.masks);
for index = 1:numel(maskNames)
    name = maskNames{index};
    observation = problem.masks.(name);
    assert(islogical(observation) && isequal(size(observation),size(problem.forcing)) && any(observation(:)), ...
        'ipm:MotionObservationMask','Observation masks must be nonempty and paired.');
    w = problem.weights(observation);
    w = w/sum(w);
    before = problem.forcing(observation)/problem.amplitudeScale;
    after = residual(observation)/problem.amplitudeScale;
    beforeL2 = sqrt(sum(w.*before.^2));
    afterL2 = sqrt(sum(w.*after.^2));
    beforeInf = max(abs(before));
    afterInf = max(abs(after));
    metrics.(name) = struct('nodeCount',nnz(observation), ...
        'originalL2',beforeL2,'modulatedL2',afterL2, ...
        'l2Ratio',afterL2/max(beforeL2,realmin), ...
        'originalInf',beforeInf,'modulatedInf',afterInf, ...
        'infRatio',afterInf/max(beforeInf,realmin));
end
end
