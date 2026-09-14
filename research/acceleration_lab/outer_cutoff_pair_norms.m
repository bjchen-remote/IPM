function report=outer_cutoff_pair_norms(first,second,radii,powers,observationNodes)
%OUTER_CUTOFF_PAIR_NORMS Compare two adaptive-grid R_X fields in one frame.
%   Both accepted states are interpolated onto the same fixed observation
%   mesh before applying the smooth outer cutoff. This measures a Cauchy
%   increment, including any intervening remesh jump, but remains subject
%   to interpolation and observation-grid refinement errors.
if nargin<3 || isempty(radii),radii=[0.05,0.1,0.2,0.4];end
if nargin<4 || isempty(powers),powers=[1,2,4,Inf];end
if nargin<5 || isempty(observationNodes),observationNodes=[401,201];end
assert(isnumeric(observationNodes) && ...
    isequal(size(observationNodes),[1,2]) && ...
    all(observationNodes>=101) && ...
    all(observationNodes==fix(observationNodes)), ...
    'ipm:OuterCutoffPairGrid','Choose two observation node counts >=101.');
for item={first,second}
    s=item{1};
    assert(all(isfield(s,{'x','y','Rx','tau'})) && ...
        isequal(size(s.Rx),[numel(s.y),numel(s.x)]) && ...
        s.x(1)<=0 && s.x(end)>=4 && s.y(1)==0 && s.y(end)>=2, ...
        'ipm:OuterCutoffPairState','Both native profiles must cover the observation box.');
end
assert(second.tau>first.tau, ...
    'ipm:OuterCutoffPairOrder','The pair must have increasing canonical time.');

xq=linspace(0,4,observationNodes(1));
yq=linspace(0,2,observationNodes(2))';
[Xq,Yq]=meshgrid(xq,yq);
a=interp2(first.x,first.y,first.Rx,Xq,Yq,'linear');
b=interp2(second.x,second.y,second.Rx,Xq,Yq,'linear');
assert(all(isfinite([a(:);b(:)])), ...
    'ipm:OuterCutoffPairInterpolation','Common-grid interpolation must be finite.');
diffField=b-a;
wx=ones(size(xq))*(xq(2)-xq(1));wx([1,end])=wx([1,end])/2;
wy=ones(size(yq))*(yq(2)-yq(1));wy([1,end])=wy([1,end])/2;
area=wy*wx;
outerX=taper(xq,3,4);
outerY=taper(yq,1,2);
outer=outerY*outerX;
distance=hypot(Yq,Xq-1);
rows=struct([]);
for radius=radii(:)'
    core=ones(size(distance));
    core(distance<=radius)=0;
    blend=distance>radius & distance<2*radius;
    core(blend)=sin(pi*(distance(blend)-radius)/(2*radius)).^2;
    mask=outer.*core;
    for domain={'bulk','wall'}
        if strcmp(domain{1},'bulk')
            weights=area;cut=mask;delta=diffField;baseline=a;latest=b;
        else
            weights=wx;cut=mask(1,:);delta=diffField(1,:);
            baseline=a(1,:);latest=b(1,:);
        end
        for p=powers(:)'
            if isinf(p)
                increment=max(abs(cut.*delta),[],'all');
                scale=max(max(abs(cut.*baseline),[],'all'), ...
                    max(abs(cut.*latest),[],'all'));
                pName='inf';
            else
                increment=sum(weights.*abs(cut.*delta).^p,'all')^(1/p);
                scale=max(sum(weights.*abs(cut.*baseline).^p,'all')^(1/p), ...
                    sum(weights.*abs(cut.*latest).^p,'all')^(1/p));
                pName=sprintf('%d',p);
            end
            row=struct('radius',radius,'domain',domain{1},'p',pName, ...
                'incrementNorm',increment, ...
                'relativeIncrementNorm',increment/max(scale,realmin));
            if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
        end
    end
end
report=struct('kind','outer_cutoff_common_grid_cauchy_increment_v1', ...
    'fromTau',first.tau,'toTau',second.tau, ...
    'observationNodes',[numel(xq),numel(yq)], ...
    'interpolation','linear_on_each_native_mesh', ...
    'rows',rows,'interpretation', ...
    ['The difference includes evolution and any remesh transfer. ' ...
     'Repeat on a finer observation grid before claiming convergence.']);
end

function values=taper(nodes,plateau,endPoint)
values=ones(size(nodes));
transition=nodes>plateau & nodes<endPoint;
values(transition)=cos(pi*(nodes(transition)-plateau)/ ...
    (2*(endPoint-plateau))).^2;
values(nodes>=endPoint)=0;
end
