function comparison = ipm_compare_bessel_two_scale_pair(referenceInput, ...
        perturbedInput,userOpts)
%IPM_COMPARE_BESSEL_TWO_SCALE_PAIR Compare matched Bessel-tail evolutions.
%   A ratio below one means the perturbed and unperturbed trajectories moved
%   closer over the tested interval; it is only evidence of attraction when
%   the decrease is resolved under longer-time and grid-refinement studies.

if nargin < 3 || isempty(userOpts)
    userOpts = struct();
end
canonicalReference = ipm.output.validate(referenceInput);
canonicalPerturbed = ipm.output.validate(perturbedInput);
reference = comparison_fields(canonicalReference);
perturbed = comparison_fields(canonicalPerturbed);
referenceMetadata = ipm_bessel_case_metadata( ...
    canonicalReference,{'centerX'});
perturbedMetadata = ipm_bessel_case_metadata( ...
    canonicalPerturbed,{'centerX'});
assert_matched_metadata(referenceMetadata,perturbedMetadata);
if ~isequal(reference.x,perturbed.x) || ~isequal(reference.y,perturbed.y)
    error('ipm:BesselPairGrid', ...
        'The matched comparison requires identical x and y grids.');
end
if abs(reference.physicalTime-perturbed.physicalTime) > ...
        100*eps(max(reference.physicalTime,perturbed.physicalTime))
    error('ipm:BesselPairTime', ...
        'The matched comparison requires identical final physical times.');
end

opts = struct('centerX',referenceMetadata.parameters.centerX, ...
    'tailXRange',[], ...
    'tailYMax',min(8,max(reference.y)));
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
centeredX = reference.x-opts.centerX;
if isempty(opts.tailXRange)
    positiveExtent = max(centeredX);
    opts.tailXRange = [max(2,0.15*positiveExtent),0.75*positiveExtent];
end
tailX = centeredX >= opts.tailXRange(1) & ...
    centeredX <= opts.tailXRange(2);
tailY = reference.y <= opts.tailYMax;

referenceInitial = first_field(canonicalReference);
perturbedInitial = first_field(canonicalPerturbed);
comparison = struct();
comparison.physicalTime = reference.physicalTime;
comparison.fullInitialDifference = weighted_relative_difference( ...
    perturbedInitial,referenceInitial,reference.x,reference.y);
comparison.fullFinalDifference = weighted_relative_difference( ...
    perturbed.physicalRho,reference.physicalRho,reference.x,reference.y);
comparison.tailInitialDifference = weighted_relative_difference( ...
    perturbedInitial(tailY,tailX),referenceInitial(tailY,tailX), ...
    reference.x(tailX),reference.y(tailY));
comparison.tailFinalDifference = weighted_relative_difference( ...
    perturbed.physicalRho(tailY,tailX), ...
    reference.physicalRho(tailY,tailX),reference.x(tailX), ...
    reference.y(tailY));
comparison.fullDifferenceRatio = comparison.fullFinalDifference / ...
    max(comparison.fullInitialDifference,eps);
comparison.tailDifferenceRatio = comparison.tailFinalDifference / ...
    max(comparison.tailInitialDifference,eps);
comparison.options = opts;
end

function fields = comparison_fields(result)
fields = struct('x',result.grid.x,'y',result.grid.y, ...
    'physicalTime',result.state.physicalTime, ...
    'physicalRho',result.physical.rho);
end

function assert_matched_metadata(reference,perturbed)
if ~strcmpi(string(reference.family),string(perturbed.family))
    error('ipm:BesselPairMetadata', ...
        'Matched Bessel runs must use the same profile family.');
end
ignored = {'perturbationAmplitude'};
referenceNames = setdiff(fieldnames(reference.parameters),ignored);
perturbedNames = setdiff(fieldnames(perturbed.parameters),ignored);
if ~isequal(sort(referenceNames),sort(perturbedNames))
    error('ipm:BesselPairMetadata', ...
        'Matched Bessel runs must expose the same base parameters.');
end
for index = 1:numel(referenceNames)
    name = referenceNames{index};
    if ~isequaln(reference.parameters.(name),perturbed.parameters.(name))
        error('ipm:BesselPairMetadata', ...
            'Matched Bessel runs differ in base parameter "%s".',name);
    end
end
end

function field = first_field(result)
if isempty(result.snapshots.rho)
    error('ipm:BesselPairSnapshots', ...
        'Matched attraction tests require storeSnapshots=true.');
end
snapshot = ipm.output.snapshotAt(result,1);
field = snapshot.physicalRho;
end

function difference = weighted_relative_difference(field,reference,x,y)
defectSquared = abs(field-reference).^2;
referenceSquared = abs(reference).^2;
difference = sqrt(trapz(y,trapz(x,defectSquared,2))) / ...
    max(sqrt(trapz(y,trapz(x,referenceSquared,2))),eps);
end
