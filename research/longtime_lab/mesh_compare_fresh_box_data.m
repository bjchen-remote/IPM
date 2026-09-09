function c=mesh_compare_fresh_box_data(a,b,window,limits)
%MESH_COMPARE_FRESH_BOX_DATA Paired physical laboratory-frame box comparison.
% Inputs are compact protocol observations; rates are already in parent units.
forward=directional(a,b,window);reverse=directional(b,a,window);c=forward;
for name=fieldnames(forward)',c.(name{1})=max(forward.(name{1}),reverse.(name{1}));end
for name={'physicalGradientMaximum','physicalQuadraticPeak','physicalCoreWidthX','physicalCoreWidthY'}
    c.([name{1} 'RelativeDifference'])=relative_inf(a.(name{1}),b.(name{1}));
end
c.cLAbsoluteDifference=abs(a.cLParentUnits-b.cLParentUnits);
c.cOmegaAbsoluteDifference=abs(a.cOmegaParentUnits-b.cOmegaParentUnits);
c.cLRelativeDifference=relative_inf(a.cLParentUnits,b.cLParentUnits);
c.cOmegaRelativeDifference=relative_inf(a.cOmegaParentUnits,b.cOmegaParentUnits);
c.passed=max(c.rhoRelativeL2,c.rhoRelativeInf)<=limits.cropRhoRelativeL2Inf && ...
    max([c.rhoXRelativeL2,c.rhoXRelativeInf,c.rhoYRelativeL2,c.rhoYRelativeInf])<=limits.cropGradientRelativeL2Inf && ...
    max(c.physicalGradientMaximumRelativeDifference,c.physicalQuadraticPeakRelativeDifference)<=limits.cropPeakRelative && ...
    max(c.physicalCoreWidthXRelativeDifference,c.physicalCoreWidthYRelativeDifference)<=limits.cropPhysicalWidthRelative && ...
    c.cLAbsoluteDifference<=max(limits.cropRateAbsolute,limits.cropRateRelative*abs(a.cLParentUnits)) && ...
    c.cOmegaAbsoluteDifference<=max(limits.cropRateAbsolute,limits.cropRateRelative*abs(a.cOmegaParentUnits));
end

function c=directional(a,b,window)
xi=find(a.x>=window(1) & a.x<=window(2));yi=find(a.y>=0 & a.y<=window(3));
assert(numel(xi)>=5 && numel(yi)>=5);[X,Y]=meshgrid(a.x(xi),a.y(yi));c=struct();
for name={'rho','rhoX','rhoY'}
    field=name{1};reference=a.(field)(yi,xi);
    sampled=interp2(b.x,b.y,b.(field),X,Y,'linear');assert(all(isfinite(sampled),'all'));
    denominator=trapz(a.y(yi),trapz(a.x(xi),reference.^2,2));assert(denominator>0);
    c.([field 'RelativeL2'])=sqrt(trapz(a.y(yi),trapz(a.x(xi),(sampled-reference).^2,2))/denominator);
    c.([field 'RelativeInf'])=relative_inf(reference,sampled);
end
end

function v=relative_inf(a,b)
v=max(abs(a-b),[],'all')/max(max(abs(a),[],'all'),realmin);
end
