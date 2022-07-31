function [PreData_temp] = PlotPredVsExp(ProfSel,OptimizationResults,numSel)
%% Extract Data
ExpData = OptimizationResults.(ProfSel).ExpData;
PreData = OptimizationResults.(ProfSel).Predata_noRes;
mz_all = OptimizationResults.(ProfSel).mz_all;
errors = OptimizationResults.(ProfSel).error;

%% Select the top number of signals to plot
% numSel: number of top peaks selected (based on the values from
% experimental profiles)

if nargin ~=3
numSel = sum(ExpData~=0);
end
if numSel>20
    numSel = 20;
end

[~,selIdx] = sort(ExpData);
selIdx = selIdx((end-numSel+1):end);
[~,ord] = sort(mz_all(selIdx));
selIdx = selIdx(ord);

ExpData_temp = ExpData(selIdx);
PreData_temp = PreData(:,selIdx);
mz_temp = mz_all(selIdx);

%% Plot experimental vs. prediction glycoprofile

% prep plotting data
plotPreData = mean(PreData_temp);
plotData = [ExpData_temp;plotPreData]';

% unbiased estimation of error (quantile)
% err = zeros(length(PreData_temp),2,2);
% PreData_temp = sort(PreData_temp,1);
% LoErrIdx = floor(size(PreData_temp,1)*0.25);
% UpErrIdx = ceil(size(PreData_temp,1)*0.75); 
% if UpErrIdx>size(PreData_temp,1); UpErrIdx = size(PreData_temp,1);end
% err(:,2,1) = plotPreData-PreData_temp(LoErrIdx,:);
% err(:,2,2) = PreData_temp(UpErrIdx,:)-plotPreData;

% biased estimation of error
err = zeros(size(PreData_temp,2),2);
err(:,2) = std(PreData_temp,[],1);
err(err ==0) = nan;

% plot bar plot
figure
hold on
rgb = {[0.25 0.25 0.25],[42 128 185]./255};
h = barwitherr(err,plotData);
xtickangle(45);

% set bar colors
legend('Experimental','Prediction','Location','northwest');
for k = 1:2
    h(k).FaceColor = rgb{k};
    h(k).EdgeColor = rgb{k};%[0 0 0];
end

% set axes, labels, and the title
xticks(1:length(mz_temp));
xticklabels(round(mz_temp));
xlabel('m/z');
ylabel('Relative Intensity');
title(sprintf(['Global Optimization for Matching Markov Model to ', ProfSel ,' Experimental Profile \n (RMSE = %0.2e)'],mean(errors)));
set(gca,'fontweight','bold')
set(gca,'TickLength',[0 0]);
hold off


end
