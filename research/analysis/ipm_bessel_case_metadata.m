function metadata = ipm_bessel_case_metadata(result,required)
%IPM_BESSEL_CASE_METADATA Read version-2 Bessel case metadata.

if nargin < 2
    required = {};
end
if isfield(result.metadata,'caseMetadata') && ...
        isstruct(result.metadata.caseMetadata) && ...
        isscalar(result.metadata.caseMetadata) && ...
        ~isempty(fieldnames(result.metadata.caseMetadata))
    metadata = result.metadata.caseMetadata;
else
    error('ipm:BesselCaseMetadata', ...
        'The version-2 result does not contain Bessel case metadata.');
end
if ~isstruct(metadata) || ~isscalar(metadata) || ...
        ~isfield(metadata,'family') || ~isfield(metadata,'parameters') || ...
        ~isstruct(metadata.parameters) || ~isscalar(metadata.parameters)
    error('ipm:BesselCaseMetadata', ...
        'The result does not contain valid Bessel case metadata.');
end
missing = required(~isfield(metadata.parameters,required));
if ~isempty(missing)
    error('ipm:BesselCaseParameters', ...
        'Missing Bessel case parameter(s): %s.',strjoin(missing,', '));
end
end
