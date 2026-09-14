function report=mesh_density_lp_oracle(distributionFile,outputFile)
%MESH_DENSITY_LP_ORACLE Fixed-node linear feasibility bound on log-width ratio.
% Research oracle only: it does not propose or transfer a solver mesh.
assert(isfile(distributionFile)&&~isfile(outputFile));
d=jsondecode(fileread(distributionFile));
x=d.x.nodes(:);y=d.y.nodes(:);
[xBound,xAxis,xDetail]=axis_bound(x,d.feature, 'x');
[yBound,yAxis,yDetail]=axis_bound(y,d.feature, 'y');
report=struct('kind','fixed_node_minimax_lp_oracle_v1', ...
    'source',distributionFile,'tau',d.tau,'nodeCount',d.nodeCount, ...
    'x',struct('baselineRatio',d.x.maximumAdjacentRatio, ...
        'oracleRatio',xBound,'axis',xAxis,'detail',xDetail), ...
    'y',struct('baselineRatio',d.y.maximumAdjacentRatio, ...
        'oracleRatio',yBound,'axis',yAxis,'detail',yDetail), ...
    'geometryOnly',true,'fullMeshCertified',false,'nativeTransferCertified',false);
fid=fopen(outputFile,'w');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end

function [bestRatio,bestAxis,detail]=axis_bound(axis,feature,direction)
n=numel(axis)-1;H=axis(end);h0=diff(axis);r0=max(exp(abs(diff(log(h0)))));
if direction=='x'
    windows=[feature.coreInterval(:)';feature.frontInterval(:)'];
    minimum=[32;20];
    anchorIndex=find(axis==1,1)-1;
    assert(~isempty(anchorIndex));
else
    windows=[0,feature.yCoreWidth];minimum=32;anchorIndex=[];
end
% Enumerate a bounded set of full-cell windows near the signed-state feature.
% This is a sufficient, slightly conservative form of the fractional count gate.
starts=cell(size(windows,1),1);
for j=1:size(windows,1)
    candidates=find(axis(1:end-minimum(j))>=windows(j,1) & ...
        axis(1+minimum(j):end)<=windows(j,2));
    if isempty(candidates)
        middle=.5*sum(windows(j,:));
        [~,q]=min(abs(axis(1:end-minimum(j))-middle));
        candidates=max(1,q-4):min(n-minimum(j)+1,q+4);
    end
    starts{j}=unique([candidates(1),candidates(round((numel(candidates)+1)/2)),candidates(end)]);
end
if direction=='y',starts{1}=1;end
combos=all_combinations(starts);
opts=optimoptions('linprog','Display','none','Algorithm','dual-simplex', ...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
bestRatio=Inf;bestAxis=[];bestRow=NaN;bestExit=NaN;
for row=1:size(combos,1)
    [A,b,Aeq,beq]=static_constraints(n,H,anchorIndex,windows,minimum,combos(row,:));
    lo=1;hi=r0+1e-6;feasible=false;bestH=[];
    for iter=1:23
        mid=.5*(lo+hi);[Atrial,btrial]=ratio_constraints(n,mid);
        [trial,~,exitflag]=linprog(zeros(n,1),[A;Atrial],[b;btrial],Aeq,beq, ...
            max(min(h0)/100,eps(H)),[],opts);
        if exitflag>0
            hi=mid;bestH=trial;feasible=true;bestExit=exitflag;
        else
            lo=mid;
        end
    end
    if feasible&&hi<bestRatio
        bestRatio=hi;bestAxis=[0;cumsum(bestH)];bestRow=row;
    end
end
if isempty(bestAxis),bestRatio=NaN;bestAxis=axis;end
detail=struct('windowStartCandidates',{starts},'combinationCount',size(combos,1), ...
    'selectedCombination',bestRow,'selectedStarts',[], ...
    'solverExitFlag',bestExit,'fullCellCounts',minimum, ...
    'intervals',windows,'anchorIndex',anchorIndex, ...
    'actualRatio',max(exp(abs(diff(log(diff(bestAxis)))))));
if isfinite(bestRow),detail.selectedStarts=combos(bestRow,:);end
end

function combos=all_combinations(starts)
if numel(starts)==1,combos=starts{1}(:);return;end
[a,b]=ndgrid(starts{1},starts{2});combos=[a(:),b(:)];
end

function [A,b,Aeq,beq]=static_constraints(n,H,anchorIndex,windows,minimum,starts)
A=zeros(2*size(windows,1),n);b=zeros(size(A,1),1);
for j=1:size(windows,1)
    i=starts(j);
    A(2*j-1,1:i-1)=1;b(2*j-1)=windows(j,2); % overwritten below
    A(2*j-1,:)=-A(2*j-1,:);b(2*j-1)=-windows(j,1);
    A(2*j,1:i+minimum(j)-1)=1;b(2*j)=windows(j,2);
end
Aeq=ones(1,n);beq=H;
if ~isempty(anchorIndex)
    row=zeros(1,n);row(1:anchorIndex)=1;
    Aeq=[Aeq;row];beq=[beq;1];
end
end

function [A,b]=ratio_constraints(n,r)
rows=(1:n-1)';A=zeros(2*(n-1),n);
for k=1:n-1
    A(k,k)=-r;A(k,k+1)=1;
    A(n-1+k,k)=1;A(n-1+k,k+1)=-r;
end
b=zeros(2*(n-1),1);
end
