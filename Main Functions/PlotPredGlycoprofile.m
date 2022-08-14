function [PreData_temp,mz_temp] = PlotPredGlycoprofile(ProfSel,OptimizationResults,numSel)
%% Extract Data
PreData = OptimizationResults.(ProfSel).Predata_noRes;
mz_all = OptimizationResults.(ProfSel).mz_all;
errors = OptimizationResults.(ProfSel).error;
leakage = mean(1-sum(PreData,2));



%% Select the top number of signals to plot
% numSel: number of top peaks selected (based on the values from
% experimental profiles)

if isempty(numSel)
    numSel = sum(ExpData~=0);
end


[~,selIdx] = sort(mean(PreData,1));
selIdx = selIdx((end-numSel+1):end);
[~,ord] = sort(mz_all(selIdx));
selIdx = selIdx(ord);

PreData_temp = PreData(:,selIdx);
ExpData_temp = zeros(1,size(PreData_temp,2));
mz_temp = mz_all(selIdx);


%% Plot experimental vs. prediction glycoprofile

% prep plotting data
plotPreData = mean(PreData_temp,1);
plotData = [ExpData_temp;plotPreData]';

% biased estimation of error
err = nan(size(PreData_temp,2),2);
err(:,2,1) = std(PreData_temp,[],1);
err(:,2,2) = std(PreData_temp,[],1);
for a = 1:length(err(:,2,2))
    if plotPreData(a)-err(a,2,1)<0
        err(a,2,1) = plotPreData(a);
    end
end

% plot bar plot
figure
hold on
rgb = {[42 128 185]./255};
h = barwitherr(err(:,2,:),plotData(:,2));
xtickangle(45);

% set bar colors
legend('Prediction','Location','northeast');
for k = 1
    h(k).FaceColor = rgb{k};
    h(k).EdgeColor = rgb{k};%[0 0 0];
end

% set axes, labels, and the title
xticks(1:length(mz_temp));
xticklabels(round(mz_temp));
xlabel('m/z');
ylabel('Relative Intensity');
title(['Predicted glycoprofile for ', strrep(ProfSel,'_','/')]);
set(gca,'fontweight','bold')
set(gca,'TickLength',[0 0]);
hold off


end