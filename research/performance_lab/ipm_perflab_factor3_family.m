function [draft,family,audit]=ipm_perflab_factor3_family(oldPolicy,oldFamily)
% Research-only five-member family. Not accepted by the native v2 contract.
assert(oldPolicy.version==2&&oldPolicy.nodeFamily.maximumTotalNodes==210000);
actual=ipm.remesh.referenceAxisFamily(oldFamily.rootX,oldFamily.rootY,oldFamily.anchor,oldPolicy);
assert(isequaln(actual,oldFamily)&&all([actual.members.resourceAdmitted]));
assert(isequal(vertcat(actual.members.cellFactors),[1,1;2,1;1,2;2,2]));
rootX=actual.rootX;rootY=actual.rootY;anchor=actual.anchor;
assert(numel(rootX)==321&&numel(rootY)==161);
xAxes={rootX,refine_x(rootX,2),refine_x(rootX,3)};
yAxes={rootY,refine_y(rootY,2)};
for k=1:4
    f=actual.members(k).cellFactors;x=xAxes{f(1)};y=yAxes{f(2)};
    assert(isequaln(x,actual.members(k).baseX)&&isequaln(y,actual.members(k).baseY)&& ...
        isequaln(ipm.mesh.quality(x,7,ipm.mesh.quadrature(x)),actual.members(k).xQuality)&& ...
        isequaln(ipm.mesh.quality(y,7,ipm.mesh.quadrature(y)),actual.members(k).yQuality));
end
x=xAxes{3};y=yAxes{2};qx=ipm.mesh.quality(x,7,ipm.mesh.quadrature(x));qy=ipm.mesh.quality(y,7,ipm.mesh.quadrature(y));
assert(hard_quality(qx)&&hard_quality(qy));
assert(isequaln(x(1:3:end),rootX)&&isequaln(y(1:2:end),rootY)&&isequal(x,-fliplr(x))&& ...
    isequal(x([1,end]),rootX([1,end]))&&isequal(y([1,end]),rootY([1,end]))&&any(x==anchor)&&any(x==-anchor));
family=actual;family.version=2;family.generator='selected_base_index_pchip_integer_v2';
family.members(5)=struct('index',5,'cellFactors',[3,2],'nodeCount',[numel(x),numel(y)], ...
    'baseX',x,'baseY',y,'resourceAdmitted',numel(x)*numel(y)<=310000,'qualityPassed',true,'xQuality',qx,'yQuality',qy);
assert(isequaln(family.members(1:4),oldFamily.members)&&isequal(family.members(5).nodeCount,[961,321]));
draft=oldPolicy;draft.version=3;draft.nodeFamily.generator=family.generator;
draft.nodeFamily.cellFactors=[1,1;2,1;1,2;2,2;3,2];draft.nodeFamily.maximumTotalNodes=310000;
draft.nodeFamily.maximumAcceptedGrowthTransitions=3;
assert(isequaln(rmfield(draft,{'version','nodeFamily'}),rmfield(oldPolicy,{'version','nodeFamily'})));
assert(isequaln(rmfield(draft.nodeFamily,{'generator','cellFactors','maximumTotalNodes','maximumAcceptedGrowthTransitions'}), ...
    rmfield(oldPolicy.nodeFamily,{'generator','cellFactors','maximumTotalNodes','maximumAcceptedGrowthTransitions'})));
counts=vertcat(family.members.nodeCount);factors=vertcat(family.members.cellFactors);
[~,order]=sortrows([prod(counts,2),(1:5)'],[1,2]);
eligible=cell(1,5);for k=1:5,eligible{k}=order(all(factors(order,:)>=factors(k,:),2));end
assert(isequal(order,[1;3;2;4;5])&&isequal(eligible{4},[4;5])&&isequal(eligible{5},5));
audit=struct('oldFourMembersEntirelyBitwise',true,'oldPolicyUnchanged',true, ...
    'draftVersion',3,'referenceFamilyVersion',2,'fullFactors',factors,'nodeCounts',counts, ...
    'nodeProducts',prod(counts,2),'ordering',order,'eligibleLevelsByCurrentLevel',{eligible}, ...
    'maximumGrowthTransitions',3,'maximumPairProposalsAllLevels',15, ...
    'oldRootKnotsExactlyPreserved',true,'allIntermediateFactorTwoKnotsPreservedByFactorThree',false, ...
    'allOriginalQualityGatesPassed',true,'oldCheckpointUpgradeable',false, ...
    'nativePolicyImplemented',false,'LUResourceQualified',false,'PDEQualification',false);
end
function x=refine_x(root,factor)
p=root(root>=0);v=pchip(linspace(0,1,numel(p)),p,linspace(0,1,factor*(numel(p)-1)+1));
v(1:factor:end)=p;v(1)=0;x=[-fliplr(v(2:end)),v];
end
function y=refine_y(root,factor)
y=pchip(linspace(0,1,numel(root)),root,linspace(0,1,factor*(numel(root)-1)+1))';y(1:factor:end)=root;y(1)=0;
end
function yes=hard_quality(q)
yes=q.quadratureWeightsStrictlyPositive&&q.maximumAdjacentCellRatio<=1.08&& ...
    q.maximumLogSpacingCurvature<=.01&&q.minimumStencilRcond>=1e-9&& ...
    q.minimumQuadratureWeightRatio>=1e-8&&q.minimumQuadratureWeightToControlWidthRatio>=.35&& ...
    q.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
