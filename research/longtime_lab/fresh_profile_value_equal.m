function same=fresh_profile_value_equal(a,b)
%FRESH_PROFILE_VALUE_EQUAL Exact stored data, including captured function data.
% Anonymous handles read from separate MAT files need not have equal handle
% identity. Compare complete FUNCTIONS descriptors, never ignore their fields.
if isa(a,'function_handle') || isa(b,'function_handle')
    same=isa(a,'function_handle') && isa(b,'function_handle');
    if same,same=fresh_profile_value_equal(functions(a),functions(b));end
elseif isstruct(a) || isstruct(b)
    same=isstruct(a) && isstruct(b) && isequal(size(a),size(b));
    if ~same,return;end
    names=sort(fieldnames(a));same=isequal(names,sort(fieldnames(b)));
    if ~same,return;end
    for k=1:numel(a)
        for n=names'
            if ~fresh_profile_value_equal(a(k).(n{1}),b(k).(n{1})),same=false;return;end
        end
    end
elseif iscell(a) || iscell(b)
    same=iscell(a) && iscell(b) && isequal(size(a),size(b));
    if ~same,return;end
    for k=1:numel(a)
        if ~fresh_profile_value_equal(a{k},b{k}),same=false;return;end
    end
else
    same=isequaln(a,b);
end
end
