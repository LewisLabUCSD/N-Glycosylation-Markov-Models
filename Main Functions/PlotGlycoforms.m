function  GlycoformData = PlotGlycoforms(ProfSel,OptimizationResults,GenericNetwork,numSel,threshold)
%% Extract Data
ExpData = OptimizationResults.(ProfSel).ExpData;
PreData = OptimizationResults.(ProfSel).Predata_raw;
Glys_raw = OptimizationResults.(ProfSel).Glys_raw;
mz_all = OptimizationResults.(ProfSel).mz_all;
AbsGlyIdx = GenericNetwork.AbsGlyIdx;
Glys = GenericNetwork.Glys;

%% Select the top number of signals to plot
% numSel: number of top peaks selected (based on the values from
% experimental profiles)

if isempty(numSel)
    numSel = 20;
    if numSel>sum(ExpData~=0)
        numSel = sum(ExpData~=0);
    end
end

[~,selIdx] = sort(ExpData);
selIdx = selIdx((end-numSel+1):end);
[~,ord] = sort(mz_all(selIdx));
selIdx = selIdx(ord);

mz_temp = mz_all(selIdx);
AbsGlyIdx = AbsGlyIdx(selIdx);
AbsGlys= cellfun(@(x) Glys(x)', AbsGlyIdx, 'UniformOutput',0);
plotData = cell(size(AbsGlys));
for a = 1:length(AbsGlys)
    glys = AbsGlys{a};
    plotData{a} = cellfun(@(x) mean(PreData(:,strcmp(Glys_raw,x))), glys);
end
plotData_norm = cellfun(@(x) x./sum(x), plotData, 'UniformOutput',0);

% Only plot glycoforms with signals > threshold
for a = 1:length(AbsGlyIdx)
    plotData_norm{a} = plotData_norm{a}(plotData{a}>threshold);
    AbsGlys{a} = AbsGlys{a}(plotData{a}>threshold);
end

keepflag = ~cellfun(@isempty,plotData_norm);
mz_temp = mz_temp(keepflag);
plotData_norm = plotData_norm(keepflag);
AbsGlys = AbsGlys(keepflag);

%% Plot experimental vs. prediction glycoprofile

% prep plot data
AllAbsGlys = [AbsGlys{:}];
plotData = zeros(length(AllAbsGlys),length(mz_temp));
for a = 1:length(plotData_norm)
    for b = 1:length(plotData_norm{a})
        plotData(strcmp(AbsGlys{a}(b),AllAbsGlys),a) = plotData_norm{a}(b);
    end
end

% Identify y labels to be highlighted
yticklabhighlight = cell(size(plotData_norm));
for a = 1:length(yticklabhighlight)
    [~,idx] = max(plotData_norm{a});
    yticklabhighlight{a} = AbsGlys{a}{idx};
end

% plot heatmap
figure;
imagesc(plotData);
colormap(flip(autumn,1));
colorbar
xlabel('m/z (Exp/Pred)','FontWeight','bold'); ylabel('Glycoforms','FontWeight','bold');
title({'Relative glycoform ratios at each m/z', '(major glycoforms highlighted red)'});
xticks(1:length(mz_temp));xticklabels(mz_temp);xtickangle(45);
yticks(1:length(AllAbsGlys));yticklabels(AllAbsGlys);

% modify axies
idxSel = ismember(AllAbsGlys,yticklabhighlight);
ticklabels = get(gca,'YTickLabel');
ticklabels_new = ticklabels;
ticklabels_new(idxSel) = strcat('\color{red} ',ticklabels_new(idxSel));
set(gca, 'YTickLabel', ticklabels_new);

% record data
GlycoformData.plotData = plotData;
GlycoformData.mz = mz_temp;
GlycoformData.AllAbsGlys = AllAbsGlys;

end