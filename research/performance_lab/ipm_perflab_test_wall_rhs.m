function report=ipm_perflab_test_wall_rhs(outputDirectory)
%IPM_PERFLAB_TEST_WALL_RHS No-LU full native 2D versus native one-line WENO.
assert(maxNumCompThreads==10&&~isfolder(outputDirectory));mkdir(outputDirectory);
registration=struct('grids',[33,17;49,25;65,33],'stretch',[0,4,10], ...
    'rates',[0,0,0,0;.23,-.17,.31,-.11;-.61,.72,-.29,.19], ...
    'boundaryModes',{{'open','closed'}},'epsilon',1e-12, ...
    'fields','Synthetic smooth density/velocities with u2(Y=0)=0. No IPM evolution or mesh admissibility claim.', ...
    'threads',10,'luCount',0,'pdeSteps',0,'originalSource',which('ipm.evolve.assembleRhs'));
records=cell(1,54);item=0;
for g=1:size(registration.grids,1)
    nx=registration.grids(g,1);ny=registration.grids(g,2);
    for stretch=registration.stretch
        xi=linspace(-1,1,nx);eta=linspace(0,1,ny)';
        if stretch==0
            x=2*xi;y=2*eta;mx=2+zeros(size(x));my=2+zeros(size(y));
        else
            x=2*sinh(stretch*xi)/sinh(stretch);mx=2*stretch*cosh(stretch*xi)/sinh(stretch);
            y=2*expm1(stretch*eta)/expm1(stretch);my=2*stretch*exp(stretch*eta)/expm1(stretch);
        end
        [X,Y]=meshgrid(x,y);
        ops=struct('nx',nx,'ny',ny,'X',X,'Y',Y,'metricX',mx,'metricY',my, ...
            'computationalSpacingX',xi(2)-xi(1),'computationalSpacingY',eta(2)-eta(1), ...
            'transportScheme','weno5_fd','wallTransportMode','advective_upwind', ...
            'spatialDiscretization','high_order','wenoEpsilon',registration.epsilon);
        Dx=ipm.mesh.fdMatrix(x,1,7);
        rho=exp(-((X-.8).^2+.5*Y.^2))+.1*sin(1.2*X).*exp(-Y);
        u1=-.8*X./(1+X.^2+Y.^2)+.15*cos(X).*(1-exp(-Y));
        u2=.2*Y./(1+X.^2+Y.^2);
        for boundary=registration.boundaryModes
            ops.transportBoundaryMode=boundary{1};
            for k=1:size(registration.rates,1)
                q=registration.rates(k,:);
                [full,base,a1]=ipm.evolve.assembleRhs(rho,u1,u2,q(1),q(2),q(3),q(4),ops,ops.transportBoundaryMode);
                [wall,wallBase,wallA1]=ipm_perflab_wall_rhs(rho(1,:),u1(1,:),q(1),q(3),q(4),ops);
                fullForcing=full*Dx';wallForcing=wall*Dx';
                fullJets=fullForcing*Dx';wallJets=wallForcing*Dx';
                item=item+1;
                records{item}=struct('nodes',[nx,ny],'stretch',stretch,'boundary',boundary{1},'rates',q, ...
                    'transportU1Exact',isequaln(a1(1,:),wallA1),'baseExact',isequaln(base(1,:),wallBase), ...
                    'rhoRateExact',isequaln(full(1,:),wall),'forcingExact',isequaln(fullForcing(1,:),wallForcing), ...
                    'forcingJetExact',isequaln(fullJets(1,:),wallJets), ...
                    'rhoRateAbsoluteInf',max(abs(full(1,:)-wall)), ...
                    'forcingAbsoluteInf',max(abs(fullForcing(1,:)-wallForcing)), ...
                    'forcingJetAbsoluteInf',max(abs(fullJets(1,:)-wallJets)), ...
                    'differentRhoRateNodes',find(full(1,:)~=wall), ...
                    'fullRhoRateWall',full(1,:),'candidateRhoRateWall',wall);
            end
        end
    end
end
assert(item==numel(records));
passed=all(cellfun(@(r)r.transportU1Exact&&r.baseExact&&r.rhoRateExact&&r.forcingExact&&r.forcingJetExact,records));
report=struct('registration',registration,'records',{records},'passed',passed, ...
    'interpretation','Original full 2D transport and candidate wall-only expression. Exact comparison includes maintained sparse first and second forcing derivatives. A failure is retained, not tolerated.');
save(fullfile(outputDirectory,'report.mat'),'report','-v7.3');
fid=fopen(fullfile(outputDirectory,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
fprintf('WALL_RHS_TINY passed=%d cases=%d exact=%d maxInf=%.6g\n',passed,item, ...
    sum(cellfun(@(r)r.rhoRateExact,records)),max(cellfun(@(r)r.rhoRateAbsoluteInf,records)));
end
