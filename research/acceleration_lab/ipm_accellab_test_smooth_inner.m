function report=ipm_accellab_test_smooth_inner(outputRoot)
%IPM_ACCELLAB_TEST_SMOOTH_INNER Polynomial/Gaussian and true time derivatives.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['smooth_inner_tiny_',token]);mkdir(destination);
p=struct('primaryTheta',.5,'holdoutTheta',[.3,.7],'power',4,'timeSteps',[.02,.01,.005],'peakAndComponentWindow',[-1,1.4], ...
    'outputDirectory',destination,'newPoissonEvaluations',0,'pdeSteps',0,'peakProjectionCount',0);
save(fullfile(destination,'registration.mat'),'p');
try
    x=linspace(-2,2,49);x=x+.012*sin(2*pi*(x+2)/4);y=linspace(0,2,33)';y=y+.01*sin(pi*y/2);
    Dx=ipm.mesh.fdMatrix(x,1,7);Dy=ipm.mesh.fdMatrix(y,1,7);window=p.peakAndComponentWindow;
    [omega,F,exact]=manufactured(x,y,0);cases=cell(1,3);thetas=[.5,.3,.7];
    for k=1:3
        theta=thetas(k);r=ipm_accellab_smooth_inner_rates(omega,F,x,y,Dx,Dy,window,theta);assert(r.valid);
        factor=sqrt(1-theta)*beta(.5,5);
        analytic=[exact.Lx*factor,exact.Ly*factor/2,exact.beta,exact.gamma,exact.velocity,exact.sigma*exact.A];
        measured=[r.rawWallWidth,r.rawVerticalWidth,r.logScaleXRate,r.logScaleYRate,r.translationRate,r.peakPrime];
        error=max(abs(measured-analytic));assert(error<1e-10);
        differences=zeros(size(p.timeSteps));
        for j=1:numel(p.timeSteps)
            h=p.timeSteps(j);[op,fp]=manufactured(x,y,h);[om,fm]=manufactured(x,y,-h);
            rp=ipm_accellab_smooth_inner_rates(op,fp,x,y,Dx,Dy,window,theta);
            rm=ipm_accellab_smooth_inner_rates(om,fm,x,y,Dx,Dy,window,theta);assert(rp.valid && rm.valid);
            finite=([rp.rawWallWidth,rp.rawVerticalWidth]-[rm.rawWallWidth,rm.rawVerticalWidth])/(2*h);
            differences(j)=max(abs(finite-[r.rawWallWidthPrime,r.rawVerticalWidthPrime]));
        end
        assert(all(differences(2:end)./differences(1:end-1)<.3));
        covariance=[];
        for scale=[.1,1,20]
            for amplitude=[.01,1,31]
                rr=ipm_accellab_smooth_inner_rates(amplitude*omega,amplitude*F,scale*x+2.3,scale*y,Dx/scale,Dy/scale,scale*window+2.3,theta);
                assert(rr.valid);
                covariance(end+1)=max(abs([rr.rawWallWidth/scale-r.rawWallWidth,rr.rawVerticalWidth/scale-r.rawVerticalWidth, ...
                    rr.logScaleXRate-r.logScaleXRate,rr.logScaleYRate-r.logScaleYRate,rr.translationRate/scale-r.translationRate])); %#ok<AGROW>
            end
        end
        assert(max(covariance)<1e-9);
        cases{k}=struct('theta',theta,'analytic',analytic,'measured',measured,'maximumAnalyticError',error, ...
            'trueParameterTimeDifferenceErrors',differences,'timeDifferenceRatios',differences(2:end)./differences(1:end-1), ...
            'translationAndUnitCovarianceErrors',covariance,'rates',r);
    end
    gaussian=cell(1,2);
    for j=1:2
        n=[49,65];n=n(j);gx=linspace(-3,3,n);gy=linspace(0,2,(n+1)/2)';
        gdX=ipm.mesh.fdMatrix(gx,1,7);gdY=ipm.mesh.fdMatrix(gy,1,7);[X,Y]=meshgrid(gx,gy);
        A=1.3;ax=.197;Lx=.65;Ly=.4;v=.07;b=-.18;g=-.23;sigma=.11;
        O=A*exp(-((X-ax)/Lx).^2-(Y/Ly).^2);OX=-2*(X-ax)/Lx^2.*O;OY=-2*Y/Ly^2.*O;
        forcing=sigma*O-(v+b*(X-ax)).*OX-g*Y.*OY;
        rows=cell(1,3);
        for k=1:3
            theta=thetas(k);r=ipm_accellab_smooth_inner_rates(O,forcing,gx,gy,gdX,gdY,window,theta);assert(r.valid);
            f=gaussian_factor(theta);expected=[Lx*f,Ly*f/2,b,g,v,sigma*A];
            measured=[r.rawWallWidth,r.rawVerticalWidth,r.logScaleXRate,r.logScaleYRate,r.translationRate,r.peakPrime];
            rows{k}=struct('theta',theta,'analytic',expected,'measured',measured,'absoluteErrors',abs(measured-expected));
        end
        gaussian{j}=struct('gridSize',[n,(n+1)/2],'results',{rows});
    end
    gaussianFineMaximum=max(cellfun(@(r)max(r.absoluteErrors),gaussian{2}.results));
    report=struct('registration',p,'polynomialCases',{cases},'gaussianCases',{gaussian}, ...
        'gaussianFineMaximumAbsoluteError',gaussianFineMaximum,'passed',gaussianFineMaximum<1e-3, ...
        'interpretation','Exact bicubic MMS, nonzero PPrime, translation/scaling covariance and finite parameter-time derivative validation. Gaussian errors are actual spatial observation errors; two grids do not certify an order. No physical IPM solve, calibration of an error gate or prior rejection change.');
    save(fullfile(destination,'report.mat'),'report');write_json(fullfile(destination,'report.json'),report);
    assert(report.passed,'ipm:SmoothGaussian','Registered Gaussian error bound failed.');
catch exception
    failure=struct('registration',p,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('SMOOTH_INNER_TINY %s pass=%d gaussianFine=%.3e\n',destination,report.passed,gaussianFineMaximum);
end
function [O,F,p]=manufactured(x,y,t)
p=struct('A',1.3*exp(.11*t),'Lx',.55*exp(-.18*t),'Ly',.4*exp(-.23*t), ...
    'a',.197+.07*t,'beta',-.18,'gamma',-.23,'velocity',.07,'sigma',.11);
[X,Y]=meshgrid(x,y);q=(X-p.a)/p.Lx;r=Y/p.Ly;
O=p.A*(1-q.^2).*(1-r.^2);OX=-2*p.A*q/p.Lx.*(1-r.^2);OY=-2*p.A*r/p.Ly.*(1-q.^2);
F=p.sigma*O-(p.velocity+p.beta*(X-p.a)).*OX-p.gamma*Y.*OY;
end
function f=gaussian_factor(theta)
z=sqrt(-log(theta));f=0;
for k=0:4
    integral=2*z;if k>0,integral=sqrt(pi/k)*erf(sqrt(k)*z);end
    f=f+nchoosek(4,k)*(-theta)^(4-k)*integral;
end
f=f/(1-theta)^4;
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
