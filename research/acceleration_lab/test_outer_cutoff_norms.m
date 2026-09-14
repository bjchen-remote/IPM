function result = test_outer_cutoff_norms()
%TEST_OUTER_CUTOFF_NORMS Check exact stationary and uniform-growth identities.
x = linspace(-5,5,101);
y = linspace(0,3,61)';
R = (1+0.2*y).*exp(-(x-2).^2);
ops = struct('x',x,'y',y,'Dx',ipm.mesh.fdMatrix(x,1,7), ...
    'remeshCount',0);
state = struct('rho',R,'ops',ops,'flow',struct('c_omega',0.13), ...
    'rhsCache',struct('rho',R,'rhoRate',zeros(size(R)), ...
        'x',x,'y',y,'remeshCount',0), ...
    'scale',struct('canonicalTime',0),'step',0);
stationary = outer_cutoff_norms(state);
rows = stationary.rows;
assert(all([rows.mixedDerivativeNorm]==0));
selected = strcmp({rows.p},'2') | strcmp({rows.p},'4');
assert(all(abs([rows(selected).candidateOuterGradientGaugeC]-0.13)<1e-11));
assert(all(abs([rows(selected).candidateGaugeRateResidual])<1e-10));

state.flow.c_omega = 0;
state.rhsCache.rhoRate = R;
growth = outer_cutoff_norms(state);
rows = growth.rows;
assert(all(abs([rows.relativeMixedDerivativeNorm]-1)<1e-11));
selected = strcmp({rows.p},'2') | strcmp({rows.p},'4');
assert(all(abs([rows(selected).candidateOuterGradientGaugeC]+1)<1e-11));
assert(all(abs([rows(selected).currentLogGradientNormRate]-1)<1e-11));
for domain = {'bulk','wall'}
    subset = rows(strcmp({rows.domain},domain{1}) & strcmp({rows.p},'2'));
    assert(all(diff([subset.gradientNorm])<=1e-11));
end
profile = struct('x',x,'y',y,'Rx',R*ops.Dx','tau',0);
same = profile;same.tau=1;
pair = outer_cutoff_pair_norms(profile,same);
assert(all([pair.rows.incrementNorm]==0));
same.Rx = 2*profile.Rx;
pair = outer_cutoff_pair_norms(profile,same);
assert(all(abs([pair.rows.relativeIncrementNorm]-0.5)<1e-10));
fine = outer_cutoff_pair_norms(profile,same,[],[],[801,401]);
assert(all(abs([fine.rows.relativeIncrementNorm]-0.5)<1e-10));
assert(isequal(fine.observationNodes,[801,401]));
result = struct('passed',true, ...
    'cases',{{'stationary','uniform_growth','common_grid_pair'}}, ...
    'radii',stationary.radii,'powers',{stationary.powers});
end
