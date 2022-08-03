function TraceGlycanSynNetwork(GenericNetwork,startGly)

fprintf('*************************************************************** \n');
fprintf('Current glycan: %s \n',startGly);

Rxns = GenericNetwork.AllrxnList_RxnTypes;
EdgeList = GenericNetwork.AllrxnList(:,[1 2]);

PreList = EdgeList(strcmp(EdgeList(:,2),startGly),1);PreList = PreList(~strcmp(startGly,PreList));
FollowList = EdgeList(strcmp(EdgeList(:,1),startGly),2);FollowList = FollowList(~strcmp(startGly,FollowList));

if ~isempty(PreList)
    fprintf('\n Preceeding glycans: \n');
    for a = 1:length(PreList)
        disp(['<a href = "matlab:TraceGlycanSynNetwork(GenericNetwork,''',sprintf('%s', PreList{a}),''')">',sprintf('%s', PreList{a}),'</a>']);
    end
end

if ~isempty(FollowList)
    fprintf('\n Following glycans: \n');
    for a = 1:length(FollowList)
        disp(['<a href = "matlab:TraceGlycanSynNetwork(GenericNetwork,''',sprintf('%s', FollowList{a}),''')">',sprintf('%s', FollowList{a}),'</a>']);
    end
end
fprintf('*************************************************************** \n');

end