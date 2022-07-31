function cmp = IdentifyBirthComp(str)

tg = {'GNb3','Ab4','GNb6)Ma6)','NNa3'};
mg = {'GNb4)Ma3','(GNb2\(\w+\)Ma6)|(GNb2Ma6)','Fa6'};

if any(cellfun(@(x) ~isempty(regexp(str,x,'once')),tg))
    cmp = 'tg';
elseif any(cellfun(@(x) ~isempty(regexp(str,x,'once')),mg))
    cmp = 'mg';
else
    if length(strfind(str,'Ma3')) == 2 && length(strfind(str,'Ma6')) == 2
        cmp = 'cg';
    else % account for ManII
        cmp = 'mg';
    end
end

end