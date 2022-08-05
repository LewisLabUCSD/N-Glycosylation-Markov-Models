function TraceGlycanSynNetwork(Prof,GenericNetwork,OptimizationResults,startGly)

fprintf('*************************************************************** \n');
fprintf('Current glycan: %s \n',startGly);

% load data
Rxns = GenericNetwork.AllrxnList_RxnTypes;
EdgeList = GenericNetwork.AllrxnList(:,[1 2]);
PseudoFlux = mean(OptimizationResults.(Prof).PseudoFlux,1);  
Rxn_flux = OptimizationResults.(Prof).Rxn_flux;  
Rxn_flux(:,1) = cellfun(@(x) x(1:end-4),Rxn_flux(:,1),'UniformOutput',false);
Rxn_flux(:,2) = cellfun(@(x) x(1:end-4),Rxn_flux(:,2),'UniformOutput',false);

% identify connecting glycans
PreList = EdgeList(strcmp(EdgeList(:,2),startGly),1);PreList = PreList(~strcmp(startGly,PreList));
FollowList = EdgeList(strcmp(EdgeList(:,1),startGly),2);FollowList = FollowList(~strcmp(startGly,FollowList));

if ~isempty(PreList)
    fprintf('\n Preceeding glycans (pseudoflux): \n');
    for a = 1:length(PreList)
        RxnSel = Rxns{strcmp(PreList{a},EdgeList(:,1)) &  strcmp(startGly,EdgeList(:,2))};
        FluxSel = PseudoFlux(strcmp(PreList{a},Rxn_flux(:,1)) &  strcmp(startGly,Rxn_flux(:,2)));
        disp(['<a href = "matlab:TraceGlycanSynNetwork(ProfSel{a},GenericNetwork,OptimizationResults,''',sprintf('%s', PreList{a}),''')">',sprintf('%s (%0.2e): %s', RxnSel, FluxSel, PreList{a}),'</a>']);
    end
end

if ~isempty(FollowList)
    fprintf('\n Following glycans (pseudoflux): \n');
    for a = 1:length(FollowList)
        RxnSel = Rxns{strcmp(FollowList{a},EdgeList(:,2)) &  strcmp(startGly,EdgeList(:,1))};
        FluxSel = PseudoFlux(strcmp(FollowList{a},Rxn_flux(:,2)) &  strcmp(startGly,Rxn_flux(:,1)));
        disp(['<a href = "matlab:TraceGlycanSynNetwork(ProfSel{a},GenericNetwork,OptimizationResults,''',sprintf('%s', FollowList{a}),''')">',sprintf('%s (%0.2e): %s', RxnSel, FluxSel,FollowList{a}),'</a>']);
    end
end
fprintf('*************************************************************** \n');

end