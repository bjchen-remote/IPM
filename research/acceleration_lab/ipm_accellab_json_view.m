function output=ipm_accellab_json_view(value)
%IPM_ACCELLAB_JSON_VIEW Non-executable JSON view; exact values remain in MAT.
if isa(value,'function_handle')
    output=struct('kind','function_handle_literal_for_json_only','literal',func2str(value), ...
        'completeExecutableHandlePreservedInMat',true);
elseif isstruct(value)
    output=value;
    for k=1:numel(value)
        for name=fieldnames(value)',output(k).(name{1})=ipm_accellab_json_view(value(k).(name{1}));end
    end
elseif iscell(value)
    output=cellfun(@ipm_accellab_json_view,value,'UniformOutput',false);
else
    output=value;
end
end
