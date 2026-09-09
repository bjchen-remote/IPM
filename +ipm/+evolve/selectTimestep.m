function [dt,stopReason,details] = selectTimestep(flow,scale,ops,time)
%IPM.EVOLVE.SELECTTIMESTEP Apply CFL and canonical/physical clock limits.
%   DETAILS records the transport rate, every candidate limit, the selected
%   limiter, and the CFL actually realized by DT.

if strcmp(ops.spatialDiscretization,'sixth_order')
    % Mapped WENO evolves on the uniform computational grid. The minimum of
    % J*dXi and the adjacent physical faces avoids optimistic CFL estimates
    % near a rapidly changing but still admitted metric.
    rate = max(abs(flow.transportU1)./ops.transportSpacingX + ...
        abs(flow.transportU2)./ops.transportSpacingY,[],'all');
else
    rate = max(abs(flow.transportU1)./ops.hx + ...
        abs(flow.transportU2)./ops.hy,[],'all');
end
if rate > 0
    dtCfl = time.cfl/rate;
else
    dtCfl = time.maxDt;
end

physicalClockSpeed = exp(scale.logC_omega-scale.logC_l);
limits = [time.maxDt,dtCfl, ...
    time.finalTime-scale.canonicalTime, ...
    (time.physicalFinalTime-scale.physicalTime)/physicalClockSpeed];
[dt,activeIndex] = min(limits);
activeLimiters = {'max_dt','cfl','canonical_final_time', ...
    'physical_final_time'};
details = struct('dt',dt,'rate',rate,'dtCfl',dtCfl, ...
    'maxDtLimit',limits(1),'canonicalLimit',limits(3), ...
    'physicalLimit',limits(4),'physicalClockSpeed',physicalClockSpeed, ...
    'realizedCfl',dt*rate, ...
    'activeLimiter',activeLimiters{activeIndex});

stopReason = '';
if dt >= time.minDt
    return;
end
if time.physicalFinalTime-scale.physicalTime <= ...
        time.minDt*physicalClockSpeed
    stopReason = 'physical_final_time';
elseif time.finalTime-scale.canonicalTime <= time.minDt
    stopReason = 'final_time';
else
    stopReason = 'time_step_below_minimum';
end
end
