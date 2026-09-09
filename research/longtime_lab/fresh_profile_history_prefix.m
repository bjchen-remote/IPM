function fresh_profile_history_prefix(old,new,replacedTerminal)
assert(isequal(sort(fieldnames(old)),sort(fieldnames(new))),'ipm:FreshCampaignPrefix','History groups changed.');
for group=fieldnames(old)'
    a=old.(group{1});b=new.(group{1});assert(isequal(sort(fieldnames(a)),sort(fieldnames(b))));
    for name=fieldnames(a)'
        x=a.(name{1});y=b.(name{1});count=numel(x)-double(replacedTerminal && ~isempty(x));
        assert(numel(y)>=count && isequaln(x(1:count),y(1:count)), ...
            'ipm:FreshCampaignPrefix','A previously accepted history prefix changed.');
        if replacedTerminal,assert(numel(x)==numel(y));end
    end
end
end
