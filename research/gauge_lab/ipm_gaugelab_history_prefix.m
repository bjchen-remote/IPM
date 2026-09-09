function [matched,audit] = ipm_gaugelab_history_prefix(parent,child)
%IPM_GAUGELAB_HISTORY_PREFIX Audit an exact grouped-history prefix.
%   Empty history groups are valid exact prefixes.  Field existence is
%   checked one scalar name at a time because ISFIELD(S,{}) does not have
%   the empty shape needed for logical indexing into FIELDNAMES(S).

groups = {'common','gauge','mesh','anisotropic'};
parentGroups = cellfun(@(name) isfield(parent,name),groups);
childGroups = cellfun(@(name) isfield(child,name),groups);
matched = isstruct(parent) && isscalar(parent) && ...
    isstruct(child) && isscalar(child) && ...
    all(parentGroups) && all(childGroups);
audit = empty_audit(matched);
if ~matched
    audit.failure = 'missing_or_nonscalar_history_group';
    return
end

for groupIndex = 1:numel(groups)
    group = groups{groupIndex};
    parentGroup = parent.(group);
    childGroup = child.(group);
    if ~isstruct(parentGroup) || ~isscalar(parentGroup) || ...
            ~isstruct(childGroup) || ~isscalar(childGroup)
        matched = false;
        audit.matched = false;
        audit.failure = 'nonscalar_history_group';
        audit.group = group;
        return
    end
    names = fieldnames(parentGroup);
    present = cellfun(@(name) isfield(childGroup,name),names);
    if ~all(present)
        matched = false;
        audit.matched = false;
        audit.failure = 'missing_history_field';
        audit.group = group;
        audit.field = strjoin(names(~present),',');
        return
    end
    for nameIndex = 1:numel(names)
        name = names{nameIndex};
        parentValues = parentGroup.(name);
        childValues = childGroup.(name);
        count = numel(parentValues);
        if numel(childValues) < count
            matched = false;
            audit.matched = false;
            audit.failure = 'short_child_history_field';
            audit = record_field(audit,group,name, ...
                parentValues,childValues);
            return
        end
        childPrefix = childValues(1:count);
        if ~isequaln(parentValues(:),childPrefix(:))
            matched = false;
            audit.matched = false;
            audit.failure = 'history_value_mismatch';
            audit = record_field(audit,group,name, ...
                parentValues,childValues);
            return
        end
    end
end
audit.matched = true;
end

function audit = empty_audit(matched)
audit = struct('matched',matched,'failure','', ...
    'group','','field','','parentCount',NaN,'childCount',NaN, ...
    'parentSize',[],'childSize',[],'parentClass','', ...
    'childClass','');
end

function audit = record_field( ...
        audit,group,name,parentValues,childValues)
audit.group = group;
audit.field = name;
audit.parentCount = numel(parentValues);
audit.childCount = numel(childValues);
audit.parentSize = size(parentValues);
audit.childSize = size(childValues);
audit.parentClass = class(parentValues);
audit.childClass = class(childValues);
end
