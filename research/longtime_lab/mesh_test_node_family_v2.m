function report=mesh_test_node_family_v2(checkpointFile,resultFile,ladderRun,oldActualRun,outDir)
%MESH_TEST_NODE_FAMILY_V2 Pure geometry/configuration migration regression.
% legacy_* are exact frozen version-one files with only the entry renamed.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
names={'ipm.config.autonomousMeshPolicy','ipm.remesh.referenceAxisFamily', ...
    'ipm.remesh.plannedAxisPairs',mfilename};
for name=names
    issues=checkcode(which(name{1}),'-id');assert(isempty(issues),jsonencode(issues));
end
files=matlab.codetools.requiredFilesAndProducts(cellfun(@which,names(1:3),'UniformOutput',false));
assert(all(contains(files,[filesep '+ipm' filesep]))&&~any(contains(files,[filesep 'research' filesep])));
report=struct('kind','registered_node_family_v2_nolu_regression','noLU',true,'noPDE',true, ...
    'sourceCheckpoint',checkpointFile,'sourceResult',resultFile,'dependencyFiles',{files}, ...
    'configurationPositive',struct([]),'configurationNegative',struct([]), ...
    'geometryPositive',struct([]),'geometryNegative',struct([]),'actualCases',struct([]));
write_json(fullfile(outDir,'registration.json'),report);
profile clear;profile on;
p1=ipm.config.autonomousMeshPolicy();
inputs={struct(),p1,struct('version',1),struct('enabled',0), ...
    struct('targetCoreCells',[42;56]),struct('targetCoreCells',[256,32]), ...
    struct('timeUnit',"NATIVE_CANONICAL")};
for k=1:numel(inputs)
    actual=ipm.config.autonomousMeshPolicy(inputs{k});expected=legacy_autonomousMeshPolicy(inputs{k});
    assert(isequaln(actual,expected)&&isequal(fieldnames(actual),fieldnames(expected)));
    report.configurationPositive=append_row(report.configurationPositive,pass(sprintf('v1_exact_%d',k)));
end
p2=ipm.config.autonomousMeshPolicy(struct('version',2,'nodeFamily',struct('maximumTotalNodes',110000)));
common=rmfield(p2,'nodeFamily');common.version=1;assert(isequaln(common,p1));
assert(isequaln(ipm.config.autonomousMeshPolicy(p2),p2));
assert(isa(p2.nodeFamily.maximumTotalNodes,'double')&&islogical(p2.nodeFamily.componentwiseNondecreasing));
report.configurationPositive=append_row(report.configurationPositive,pass('v2_common_exact_and_roundtrip'));
canonical=p2;canonical.nodeFamily.generator=upper(string(canonical.nodeFamily.generator));
canonical.nodeFamily.ordering=upper(string(canonical.nodeFamily.ordering));
canonical.nodeFamily.componentwiseNondecreasing=1;
canonical.nodeFamily.maximumTotalNodes=uint32(110000);
assert(isequaln(ipm.config.autonomousMeshPolicy(canonical),p2));
report.configurationPositive=append_row(report.configurationPositive,pass('v2_normalizes_registered_input_types'));
choices=struct('rescalingMode','dynamic','dynamicScaleGeometry','isotropic', ...
    'symmetryMode','double_odd_omega','lengthGauge','transport_anchor', ...
    'cOmegaGauge','wall_omega_quadratic_peak','spatialDiscretization','high_order', ...
    'transportScheme','weno5_fd','timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'adaptiveRemesh',true,'initialAnalyticRemesh',true,'adaptiveLevels',[.1,.5,.9]);
assert(isequaln(ipm.config.autonomousMeshPolicy(p2,choices),p2));
disabled=p2;disabled.enabled=false;
assert(isequaln(ipm.config.autonomousMeshPolicy(disabled,struct()),disabled));
report.configurationPositive=append_row(report.configurationPositive,pass('v2_enabled_and_disabled_combination_contract'));
incompatible=choices;incompatible.initialAnalyticRemesh=false;
report.configurationNegative=append_row(report.configurationNegative,reject('v2_missing_analytic_initialization', ...
    @()ipm.config.autonomousMeshPolicy(p2,incompatible),'ipm:AutonomousMeshCombination'));
incompatible=choices;incompatible.dynamicScaleGeometry='anisotropic';
report.configurationNegative=append_row(report.configurationNegative,reject('v2_incompatible_geometry', ...
    @()ipm.config.autonomousMeshPolicy(p2,incompatible),'ipm:AutonomousMeshCombination'));
badInputs={struct('version',2),struct('version',2,'nodeFamily',struct()), ...
    struct('nodeFamily',struct('maximumTotalNodes',110000)), ...
    struct('version',0),struct('version',3),struct('version',true)};
badLabels={'v2_missing_family','v2_missing_cap','v1_family_forbidden','version_zero','version_three','logical_version'};
for value={0,-1,.5,Inf,NaN,[110000,110000],"110000",true,complex(110000,1)}
    bad=p2;bad.nodeFamily.maximumTotalNodes=value{1};badInputs{end+1}=bad; %#ok<AGROW>
    badLabels{end+1}=sprintf('invalid_cap_%d',numel(badLabels)-5); %#ok<AGROW>
end
for pair={{'generator','other'},{'ordering','registration_index'}, ...
        {'cellFactors',[1,1;1,2;2,1;2,2]},{'cellFactors',[1,1,2,1,1,2,2,2]}, ...
        {'componentwiseNondecreasing',false},{'maximumAcceptedGrowthTransitions',3},{'unknown',1}}
    bad=p2;bad.nodeFamily.(pair{1}{1})=pair{1}{2};badInputs{end+1}=bad; %#ok<AGROW>
    badLabels{end+1}=['fixed_family_',pair{1}{1}]; %#ok<AGROW>
end
bad=p2;bad.qualityLimits.maxAdjacentCellRatio=1.09;badInputs{end+1}=bad;badLabels{end+1}='quality_not_relaxed';
bad=p2;bad.search.fineCells=[32,64];badInputs{end+1}=bad;badLabels{end+1}='search_not_changed';
for k=1:numel(badInputs)
    value=badInputs{k};report.configurationNegative=append_row(report.configurationNegative, ...
        reject(badLabels{k},@()ipm.config.autonomousMeshPolicy(value),'ipm:BadAutonomousMeshPolicy'));
end

cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;r=ipm.output.validate(resultFile);
assert(isequaln(s.rho,r.state.rho)&&isequal(s.x,r.grid.x)&&isequal(s.y,r.grid.y)&& ...
    s.step==r.state.steps&&s.scale.canonicalTime==r.state.canonicalTime&& ...
    s.scale.physicalTime==r.state.physicalTime&&s.normalizedTime==r.state.normalizedTime&& ...
    isequaln(cp.payload.log.history,r.history)&&isequaln(s.config,r.config)&& ...
    s.remeshCount==r.grid.remeshCount&&strcmp(s.runMetadata.caseId,r.metadata.caseId)&& ...
    r.scale.Cx==exp(s.scale.logC_l)&&r.scale.Cy==exp(s.scale.logC_l)&& ...
    r.scale.Comega==exp(s.scale.logC_omega)&&r.scale.Xshift==s.scale.X_shift&&all(ipm.output.trustedMask(r)));
assert(s.step==1134&&isequal([numel(s.x),numel(s.y)],[321,161]));
rootX=s.baseX;rootY=s.baseY;anchor=s.config.scaling.transportAnchorX;
assert(isequal(rootX,s.runMetadata.autonomousMesh.initialization.selectedBaseX)&& ...
    isequal(rootY,s.runMetadata.autonomousMesh.initialization.selectedBaseY));
family=ipm.remesh.referenceAxisFamily(rootX,rootY,anchor,p2);
assert(isequal([family.members.index],1:4)&& ...
    isequal(vertcat(family.members.nodeCount),[321,161;641,161;321,321;641,321])&& ...
    isequal([family.members.resourceAdmitted],[true,true,true,false])&&all([family.members.qualityPassed]));
assert(isequal(family.rootX,rootX)&&isequal(family.rootY,rootY));
counts=prod(vertcat(family.members.nodeCount),2);
[~,order]=sortrows([counts,[family.members.index]']);assert(isequal(order(:)',[1,3,2,4]));
for k=1:4
    m=family.members(k);f=m.cellFactors;
    assert(isequal(m.baseX(1:f(1):end),rootX)&&isequal(m.baseY(1:f(2):end),rootY)&& ...
        any(m.baseX==anchor)&&any(m.baseX==-anchor)&&isequal(m.baseX,-fliplr(m.baseX))&& ...
        isequal(m.baseX([1,end]),rootX([1,end]))&&isequal(m.baseY([1,end]),rootY([1,end])));
    assert(isequaln(m.xQuality,ipm.mesh.quality(m.baseX,7,ipm.mesh.quadrature(m.baseX)))&& ...
        isequaln(m.yQuality,ipm.mesh.quality(m.baseY,7,ipm.mesh.quadrature(m.baseY))));
end
report.geometryPositive=append_row(report.geometryPositive,pass('all_registered_members_quality_and_exact_knots'));
assert(isequaln(family,ipm.remesh.referenceAxisFamily(rootX,rootY,anchor,p2)));
report.geometryPositive=append_row(report.geometryPositive,pass('family_deterministic_full_reconstruction'));
for factor=[.37,4]
    units=ipm.remesh.referenceAxisFamily(factor*rootX,factor*rootY,factor*anchor,p2);relative=0;
    for k=1:4
        relative=max(relative,max([max(abs(units.members(k).baseX/factor-family.members(k).baseX))/max(abs(rootX)), ...
            max(abs(units.members(k).baseY/factor-family.members(k).baseY))/max(abs(rootY))]));
    end
    assert(relative<1e-14);
    report.geometryPositive=append_row(report.geometryPositive,pass(sprintf('axis_unit_covariance_%g',factor)));
end
rootOnly=p2;rootOnly.nodeFamily.maximumTotalNodes=321*161;
rootFamily=ipm.remesh.referenceAxisFamily(rootX,rootY,anchor,rootOnly);
assert(isequal([rootFamily.members.resourceAdmitted],[true,false,false,false])&& ...
    all([rootFamily.members.qualityPassed]));
report.geometryPositive=append_row(report.geometryPositive,pass('unadmitted_members_still_checked'));
bad=p2;bad.nodeFamily.maximumTotalNodes=321*161-1;
report.geometryNegative=append_row(report.geometryNegative,reject('root_over_resource_cap', ...
    @()ipm.remesh.referenceAxisFamily(rootX,rootY,anchor,bad),'ipm:AutonomousMeshResourceCap'));
report.geometryNegative=append_row(report.geometryNegative,reject('v1_family_forbidden', ...
    @()ipm.remesh.referenceAxisFamily(rootX,rootY,anchor,p1),'ipm:AutonomousMeshReferenceFamily'));
report.geometryNegative=append_row(report.geometryNegative,reject('inexact_anchor', ...
    @()ipm.remesh.referenceAxisFamily(rootX,rootY,anchor+eps(anchor),p2),'ipm:AutonomousMeshReferenceFamily'));
badX=rootX;badX(end)=badX(end)+eps(badX(end));
report.geometryNegative=append_row(report.geometryNegative,reject('broken_symmetry', ...
    @()ipm.remesh.referenceAxisFamily(badX,rootY,anchor,p2),'ipm:AutonomousMeshReferenceFamily'));
badX=rootX;positive=find(badX>anchor,1);badX(positive)=badX(positive-1)+.001*(badX(positive)-badX(positive-1));
badX(numel(badX)+1-positive)=-badX(positive);
report.geometryNegative=append_row(report.geometryNegative,reject('axis_quality_even_when_growth_unadmitted', ...
    @()ipm.remesh.referenceAxisFamily(badX,rootY,anchor,rootOnly),'ipm:AutonomousMeshReferenceQuality'));

view=accepted_view(struct('rho',s.rho,'x',s.x,'y',s.y));
reference=struct('x',rootX,'y',rootY);
[legacy,legacyReport]=legacy_plannedAxisPairs(view,reference,anchor,p1);
[current,currentReport]=ipm.remesh.plannedAxisPairs(view,reference,anchor,p1);
assert(isequaln(legacy,current)&&isequaln(strip_stacks(legacyReport),strip_stacks(currentReport)));
report.actualCases=append_row(report.actualCases,case_row('strict_1134_v1_bitwise',currentReport,true));
report.sourceNativeStrict=true;report.sourcePairingExact=true;report.sourceStep=s.step;
% The old experiment was in node-product order; production indices are fixed.
oldLevels=[1,3,2,4];wide=p2;wide.nodeFamily.maximumTotalNodes=250000;
for index=1:4
    member=family.members(index);reference=struct('x',member.baseX,'y',member.baseY);
    [candidate,axisReport]=ipm.remesh.plannedAxisPairs(view,reference,anchor,wide);
    old=load(fullfile(ladderRun,sprintf('level_%d_axes.mat',oldLevels(index))));
    exact=isequaln(candidate,old.axes)&&isequaln(strip_stacks(axisReport.xTrials),strip_stacks(old.axisReport.xTrials))&& ...
        isequaln(strip_stacks(axisReport.yTrials),strip_stacks(old.axisReport.yTrials))&& ...
        isequal(reference,old.targetReference);
    assert(exact&&numel(candidate)<=3&&axisReport.sameBoxSameNodeCount==(index==1));
    report.actualCases=append_row(report.actualCases,case_row(sprintf('strict_1134_family_%d',index),axisReport,exact));
    save(fullfile(outDir,sprintf('member_%d_axes.mat',index)),'candidate','axisReport','member','-v7.3');
    fprintf('NODE_FAMILY_ACTUAL index=%d N=%dx%d x=%d y=%d pair=%d exact=1\n',index,member.nodeCount,nnz([axisReport.xTrials.admissible]),nnz([axisReport.yTrials.admissible]),numel(candidate));
end
over=family.members(4);reference=struct('x',over.baseX,'y',over.baseY);
report.geometryNegative=append_row(report.geometryNegative,reject('planner_target_over_cap', ...
    @()ipm.remesh.plannedAxisPairs(view,reference,anchor,p2),'ipm:AutonomousMeshResourceCap'));
reference=struct('x',family.members(2).baseX,'y',rootY);
report.geometryNegative=append_row(report.geometryNegative,reject('v1_differentN_rejected', ...
    @()ipm.remesh.plannedAxisPairs(view,reference,anchor,p1),''));
for label={'original_t0','full_box_step3326','fresh_step3668'}
    old=load(fullfile(oldActualRun,[label{1},'.mat']));view=accepted_view(old.snapshot);
    [legacy,legacyReport]=legacy_plannedAxisPairs(view,old.reference,old.anchor,p1);
    [current,currentReport]=ipm.remesh.plannedAxisPairs(view,old.reference,old.anchor,p1);
    exact=isequaln(legacy,current)&&isequaln(strip_stacks(legacyReport),strip_stacks(currentReport));assert(exact);
    report.actualCases=append_row(report.actualCases,case_row([label{1},'_v1_bitwise'],currentReport,exact));
end
% A broad analytic fixture tests that a valid same-N keep is never returned
% in place of a requested larger reference. It is not native PDE evidence.
x=linspace(-4,4,321);y=linspace(0,1,129)';[X,Y]=meshgrid(x,y);
view=accepted_view(struct('rho',(erf((X-2)/1.8)-erf((X+2)/1.8)).*exp(-Y.^2),'x',x,'y',y));
tinyFamily=ipm.remesh.referenceAxisFamily(x,y,1,wide);
[keep,~]=ipm.remesh.plannedAxisPairs(view,struct('x',x,'y',y),1,p1);assert(~isempty(keep)&&keep(1).unchanged);
[grown,grownReport]=ipm.remesh.plannedAxisPairs(view,struct('x',tinyFamily.members(2).baseX,'y',y),1,wide);
assert(~any([grown.unchanged])&&all(arrayfun(@(v)numel(v.x)==641,grown))&& ...
    isequal(grownReport.nodeCount,[641,129])&&~grownReport.sameBoxSameNodeCount);
report.geometryPositive=append_row(report.geometryPositive,pass('sameN_keep_excluded_from_differentN_target'));
report.broadFixtureGrowthStatus=grownReport.status;
profile off;profileInfo=profile('info');profileNames={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.evolve.advance','ipm.solve'}
    assert(~any(strcmp(profileNames,forbidden{1})),'ipm:NodeFamilyUnexpectedPDE','No solver or large LU work is permitted.');
end
report.allPassed=all([report.configurationPositive.passed])&&all([report.configurationNegative.passed])&& ...
    all([report.geometryPositive.passed])&&all([report.geometryNegative.passed])&&all([report.actualCases.passed]);
report.nodeProductOrder=order(:)';report.productionDependencyCount=numel(files);
save(fullfile(outDir,'report.mat'),'report','family','profileInfo','-v7.3');write_json(fullfile(outDir,'report.json'),report);
assert(report.allPassed);
fprintf('NODE_FAMILY_V2_COMPLETE configPositive=%d configNegative=%d geometryPositive=%d geometryNegative=%d actual=%d noLU=1\n', ...
    numel(report.configurationPositive),numel(report.configurationNegative),numel(report.geometryPositive),numel(report.geometryNegative),numel(report.actualCases));
end

function view=accepted_view(s)
Dx=ipm.mesh.fdMatrix(s.x,1,7);view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',s.rho*Dx','trusted',true);
end
function value=strip_stacks(value)
if ~isstruct(value),return;end
for k=1:numel(value)
    for field=fieldnames(value)'
        if strcmp(field{1},'failure')&&isstruct(value(k).failure)&&isfield(value(k).failure,'stack')
            value(k).failure=rmfield(value(k).failure,'stack');
        elseif isstruct(value(k).(field{1}))
            value(k).(field{1})=strip_stacks(value(k).(field{1}));
        end
    end
end
end
function row=case_row(label,r,passed)
row=struct('label',label,'nodeCount',r.nodeCount,'admittedX',nnz([r.xTrials.admissible]), ...
    'admittedY',nnz([r.yTrials.admissible]),'candidateCount',numel(r.selectedPairs),'passed',passed);
end
function row=pass(label)
row=struct('label',label,'passed',true);
end
function row=reject(label,f,expected)
actual='';message='';
try
    f();
catch e
    actual=e.identifier;message=e.message;
end
passed=~isempty(actual)&&(isempty(expected)||strcmp(actual,expected));
row=struct('label',label,'expected',expected,'actual',actual,'message',message,'passed',passed);assert(passed,label);
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
