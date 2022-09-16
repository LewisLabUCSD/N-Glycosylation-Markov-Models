function [FluxesbyComp,Subnetwork] = PlotModelFluxes(ProfSel,OptimizationResults,GenericNetwork,numSel)
%% Extract & Prep Data 
RxnTypes = GenericNetwork.RxnTypes;
RxnTypes_noBr = cellfun(@(x) regexprep(x,'\_B\d',''),RxnTypes,'UniformOutput',0);
AllrxnList = GenericNetwork.AllrxnList;
AllrxnList_RxnTypes = GenericNetwork.AllrxnList_RxnTypes;
PseudoFlux = OptimizationResults.(ProfSel).PseudoFlux;
Rxn_flux = OptimizationResults.(ProfSel).Rxn_flux;
PseudoConc = OptimizationResults.(ProfSel).PseudoConc;
Glys_conc = OptimizationResults.(ProfSel).Glys_conc;
Glys = GenericNetwork.Glys;

% Identify RxnType compartments
AllComp = cellfun(@(x) AllrxnList{find(strcmp(AllrxnList(:,3),x),1),4}(1:4), RxnTypes_noBr, 'UniformOutput',0);
Comp = unique(AllComp);
sets = cellfun(@(x) find(strcmp(x,AllComp)), Comp, 'UniformOutput',0);
setnames = {'cis','med','trans'};

%% Visualize pseudo-fluxes by reaction types

% compute pseudo-fluxes by reaction types
FluxesbyComp = zeros(size(PseudoFlux,1),length(RxnTypes));
for a = 1:size(FluxesbyComp,1)
    for b = 1:size(FluxesbyComp,2)
        FluxesbyComp(a,b) = sum(PseudoFlux(a,strcmp(AllrxnList_RxnTypes,RxnTypes{b})));
    end
end
FluxesbyComp_mean = mean(FluxesbyComp,1);

% biased estimation of error
err = zeros(1,size(FluxesbyComp,2),2);
err(1,:,1) = std(FluxesbyComp,[],1);
err(1,:,2) = std(FluxesbyComp,[],1);
for a = 1:length(err(1,:,2))
    if FluxesbyComp_mean(a)-err(1,a,1)<0
        err(1,a,1) = FluxesbyComp_mean(a);
    end
end

% Plot pseudo-fluxes by compartments
figure
hold on
for k = 1:3
    
    subplot(1,3,k)
    hold on
    h = barwitherr(err(:,sets{k},:), FluxesbyComp_mean(sets{k}));
    xlabel(['Rxn Types (',setnames{k},')']);
    ylabel('Total relative pseudo-fluxes');
    xticks(1:length(RxnTypes(sets{k})));
    xticklabels(strrep(RxnTypes(sets{k}),'_',' '));
    xtickangle(45);
    set(gca,'fontweight','bold')
    set(gca,'TickLength',[0 0]);

    % set bar colors
    h(1).FaceColor = [42 128 185]./255;
    hold off
    
end
hold off
sgtitle(['Model pseudo-fluxes through reaction types (',ProfSel,')']);
hold off

%% Visualize reactant pseudo-concentrations by reaction types

% compute pseudo-fluxes by reaction types
RctPseudoConcbyComp = zeros(size(PseudoConc,1),length(RxnTypes));
for a = 1:size(RctPseudoConcbyComp,1)
    for b = 1:size(RctPseudoConcbyComp,2)
        RctSel_idx = find(strcmp(AllrxnList_RxnTypes,RxnTypes{b}));
        RctSel = strcat(AllrxnList(RctSel_idx,1),cellfun(@(x) x(1:4), AllrxnList(RctSel_idx,4),'UniformOutput',0));
        RctPseudoConcbyComp(a,b) = sum(cellfun(@(x) PseudoConc(a,strcmp(x,Glys_conc)), RctSel));
    end
end
RctPseudoConcbyComp_mean = mean(RctPseudoConcbyComp,1);

% biased estimation of error
err = zeros(1,size(RctPseudoConcbyComp,2),2);
err(1,:,1) = std(RctPseudoConcbyComp,[],1);
err(1,:,2) = std(RctPseudoConcbyComp,[],1);
for a = 1:length(err(1,:,2))
    if RctPseudoConcbyComp_mean(a)-err(1,a,1)<0
        err(1,a,1) = RctPseudoConcbyComp_mean(a);
    end
end

% Plot pseudo-fluxes by compartments
figure
hold on
for k = 1:3
    
    subplot(1,3,k)
    hold on
    h = barwitherr(err(:,sets{k},:), RctPseudoConcbyComp_mean(sets{k}));
    xlabel(['Rxn Types (',setnames{k},')']);
    ylabel('Total relative pseudo-fluxes');
    xticks(1:length(RxnTypes(sets{k})));
    xticklabels(strrep(RxnTypes(sets{k}),'_',' '));
    xtickangle(45);
    set(gca,'fontweight','bold')
    set(gca,'TickLength',[0 0]);

    % set bar colors
    h(1).FaceColor = [42 128 185]./255;
    hold off
    
end
hold off
sgtitle(['Model reactant pseudo-concentrations for reaction types (',ProfSel,')']);
hold off


%% Visualize pseudo-fluxes in network

%  Self-loop edges (absorption rxns) ignored
SelIdx = ~contains(Rxn_flux(:,1),'[ab]');
Rxn_flux = Rxn_flux(SelIdx,:);PseudoFlux = PseudoFlux(:,SelIdx);
PseudoFlux_mean = mean(PseudoFlux,1);
SelIdx = contains(Rxn_flux(:,1),'[tg]') & contains(Rxn_flux(:,2),'[ab]');
AbsGlySel = Rxn_flux(SelIdx,2);
PseudoFlux_mean_AbsGlySel = PseudoFlux_mean(SelIdx);

%%%%%%%%%%%%%%%%%%%%%%%% Plot background network %%%%%%%%%%%%%%%%%%%%%%%%
G = digraph(Rxn_flux(:,1),Rxn_flux(:,2));
figure('Units','normalized','Position',[0 0 1 1]);
hold on
sgtitle({['Major Model Fluxes for ',ProfSel],'(major edges highlited blue and absoprtion glycans highlighted red)'},'FontWeight','bold');
subplot(1,2,1);
hold on
h = plot(G,'LayOut','layered');
h.EdgeColor = '#ecf0f1';
h.NodeColor = '#bec3c7';
h.NodeLabel = {};

%%%%%%%%%%%%%%%%%%%%%%%% Plot trimmed network %%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~=4
   numSel = 20;
end

[~,idx] = sort(PseudoFlux_mean_AbsGlySel,'descend');
AbsGlySel = AbsGlySel(idx(1:numSel));
[SubNodes, SubEdges] = extractIntermediateNetwork(G, AbsGlySel);
SubNodes = setdiff(SubNodes,AbsGlySel);

% highlihgt on the figure
highlight(h,AbsGlySel,'NodeColor','#D95319','MarkerSize',4);
highlight(h,SubEdges(:,1), SubEdges(:,2), 'EdgeColor',[42 128 185]./255);
highlight(h,SubNodes,'NodeColor',[42 128 185]./255,'MarkerSize',4);
title('Complete Network');
set(gca,'xtick',[]);set(gca,'ytick',[]);
hold off

% Construct graph
G_sub = digraph(SubEdges(:,1),SubEdges(:,2));

% compute edge weights
EdgeWt = zeros(size(G_sub.Edges.EndNodes,1),1);
for a = 1:length(EdgeWt)
    EdgeWt(a) = PseudoFlux_mean(strcmp(Rxn_flux(:,1),G_sub.Edges.EndNodes(a,1)) & strcmp(Rxn_flux(:,2),G_sub.Edges.EndNodes(a,2)));
end

% Compute NodeWeights 
NodeWt = cellfun(@(x) mean(PseudoConc(:,strcmp(x,Glys_conc))),G_sub.Nodes.Name);
NodeWt = normalize(NodeWt,'range',[2 8]);

AllPaths = cell(1,length(AbsGlySel));
MajorEdges_temp = {};
for a = 1:length(AllPaths)
    AllPaths{a} = ConcentrationGradientChaser(G_sub,NodeWt,AbsGlySel{a});
    MajorEdges_temp = [MajorEdges_temp;strcat(AllPaths{a}(1:end-1)','->',AllPaths{a}(2:end)')];
end
MajorEdges_temp = unique(MajorEdges_temp);
MajorEdges = cell(length(MajorEdges_temp),2);
for a = 1:length(MajorEdges_temp)
    edge = strsplit(MajorEdges_temp{a},'->');
    MajorEdges(a,:) = edge;
end

% Plot subtracted network
subplot(1,2,2);
hold on
h = plot(G_sub,'LayOut','layered','EdgeColor',[42 128 185]./255,'NodeColor',[42 128 185]./255);
highlight(h,AbsGlySel,'NodeColor','#D95319');
for a = 1:length(AllPaths)
    highlight(h,AllPaths{a},'EdgeColor','r','LineWidth',2)
end
h.MarkerSize = NodeWt;
title('Trimmed Network');
h.NodeLabel = {};
set(gca,'xtick',[]);set(gca,'ytick',[]);
hold off

hold off

% Record data
Subnetwork.G_trimmed = G_sub;
Subnetwork.NodeWt_sub = NodeWt;
Subnetwork.EdgeWt_sub = EdgeWt;
Subnetwork.RedEdges_sub = MajorEdges;

% %%%%%%%%%%%%%%% Plot major branching points %%%%%%%%%%%%

% Get all edges from AllPaths
MajorEdges = {};
for a = 1:length(AllPaths)
    MajorEdges = [MajorEdges;[AllPaths{a}(1:end-1)',AllPaths{a}(2:end)']];
end
MajorEdges = cellstr(unique(categorical(MajorEdges), 'rows'));
NodeNames = unique([MajorEdges(:,1);MajorEdges(:,2)]);

% Remove points with bot in and out degrees = 1
MajorEdges = RemoveNonBranchingNodes(MajorEdges);

% Get edge names (RxnTypes)
[row,col] = ind2sub(size(GenericNetwork.TM), GenericNetwork.AllrxnList_TMidx);
rowGly = Glys(row);colGly = Glys(col);
EdgeNames = cell(size(MajorEdges,1),1);
EdgeFluxes = zeros(size(MajorEdges,1),1);
for a = 1:length(EdgeNames)
    idx = find(strcmp(rowGly,MajorEdges{a,1}) & strcmp(colGly,MajorEdges{a,2}),1);
    if ~isempty(idx)
    EdgeNames{a} = AllrxnList_RxnTypes{idx};
    EdgeFluxes(a) = mean(PseudoFlux(:,idx));
    else % used outNode psuedo-concentration
        EdgeNames{a} = '';
        EdgeFluxes(a) = PseudoConc(strcmp(MajorEdges{a,2},Glys_conc));
    end
end

% Construct digraph
G = digraph(MajorEdges(:,1),MajorEdges(:,2),EdgeFluxes);

% Get Node weights
nodes = G.Nodes.Name;
NodeWt = cellfun(@(x) mean(PseudoConc(:,strcmp(x,Glys_conc))),nodes);
NodeWt = normalize(log2(NodeWt),'range',[2,15]);

% Plot 
figure
hold on
h = plot(G,'EdgeColor','r','LineWidth',2);
h.NodeLabel = {};
h.MarkerSize = NodeWt;
labeledge(h,MajorEdges(:,1),MajorEdges(:,2),strrep(EdgeNames,'_',' '));
highlight(h,AbsGlySel,'NodeColor','#D95319');
xticks([]);yticks([]);
title(['Major intermediate glycans in N-glycosylation biosynthetic network', '(',ProfSel,')']);
% title({'Major branching points in N-glycosylation biosynthetic network',[ProfSel, ' ,Node size proportional to log2(Pseudo-Concentration)']});
hold off

% Store varaibles
OptimizationResults.ProfSel.Subnetwork.G_branching = G;
OptimizationResults.ProfSel.Subnetwork.NodeWt_branching = NodeWt;

end