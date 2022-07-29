function [data,xval_tot]= PlotTPbyComp(ProfSel,OptimizationResults,GenericNetwork)

% load data
xval_tot = OptimizationResults.(ProfSel).xval(~isoutlier(OptimizationResults.(ProfSel).fval),:);
RxnTypes = GenericNetwork.RxnTypes;
RxnTypes_noBr = cellfun(@(x) regexprep(x,'\_B\d',''),RxnTypes,'UniformOutput',0);
AllrxnList = GenericNetwork.AllrxnList;
StericFlag = OptimizationResults.(ProfSel).OptimizationProblem.StericFlag;
AppliedGeneidx = OptimizationResults.(ProfSel).OptimizationProblem.AppliedGeneidx;
stericRxns = OptimizationResults.(ProfSel).OptimizationProblem.stericRxns;

% Identify RxnType compartments
AllComp = cellfun(@(x) AllrxnList{find(strcmp(AllrxnList(:,3),x),1),4}(1:4), RxnTypes_noBr, 'UniformOutput',0);
Comp = unique(AllComp);
setnames = {'cis','med','trans'};
[~,transportIdx] = ismember({'cg2mg','mg2tg','tg2ab'},RxnTypes);

if StericFlag
    AllComp = [AllComp;repmat({'[SF]'},length(stericRxns),1)];
    setnames = [setnames,{'SF'}];
    Comp = [Comp;{'[SF]'}];
end
sets = cellfun(@(x) find(strcmp(x,AllComp)), Comp, 'UniformOutput',0);

% Normalize scaling factors by compartments
xval_tot = 10.^xval_tot;
for k1 = 1:size(xval_tot,1)
    for k2 = 1:3
        vec = xval_tot(k1,:);
        vec(sets{k2}) = vec(sets{k2})./vec(transportIdx(k2));
    end
    xval_tot(k1,:) = vec;
end

% Calculate variables
xval_tot = log10(xval_tot);
data = median(xval_tot,1);


%% Plot data
figure
hold on
for k = 1:length(Comp)
    
    subplot(1,length(Comp),k)
    hold on
    boxplot(xval_tot(:,sets{k}),'Colors',[42 128 185]./255,'Symbol','ro','OutlierSize',2,'MedianStyle','target');
    
    ylabel('Log10 Relative Transition Probabilities');
    if strcmp(Comp{k},'[SF]')
        xticks(1:length(stericRxns));
        xticklabels(stericRxns);
        xlabel('Steric Factors');
    else
        xticks(1:length(RxnTypes(sets{k})));
        xticklabels(strrep(RxnTypes(sets{k}),'_',' '));
        xlabel(['Rxn Types (',setnames{k},')']);
    end
    xtickangle(45);
    set(gca,'fontweight','bold')
    set(gca,'TickLength',[0 0]);
    hold off
    
end
hold off

sgtitle(['Transition Probabilities by Compartment (',ProfSel,')']);
hold off



end