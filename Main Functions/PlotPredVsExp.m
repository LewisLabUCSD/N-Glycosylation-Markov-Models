function [PreData_temp,mz_temp] = PlotPredVsExp(ProfSel,OptimizationResults,numSel)
%% Extract Data
ExpData = OptimizationResults.(ProfSel).ExpData;
PreData = OptimizationResults.(ProfSel).Predata_noRes;
mz_all = OptimizationResults.(ProfSel).mz_all;
errors = OptimizationResults.(ProfSel).error;
leakage = mean(1-sum(PreData,2));
randomModelflag = ismember({'RandomResults'},fieldnames(OptimizationResults.(ProfSel)));

%% Select the top number of signals to plot
% numSel: number of top peaks selected (based on the values from
% experimental profiles)

if isempty(numSel)
    numSel = sum(ExpData~=0);
end
if numSel>sum(ExpData~=0)
    numSel = sum(ExpData~=0);
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
plotPreData = mean(PreData_temp,1);
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
rgb = {[0.25 0.25 0.25],[42 128 185]./255};

if all(size(plotPreData) == [1,1])
    h = bar(1,plotData(1));h2 = bar(2,plotData(2));

    h.FaceColor = rgb{1};h.EdgeColor = rgb{1};
    h2.FaceColor = rgb{2};h2.EdgeColor = rgb{2};
    errorbar(2,plotPreData,err(1,2,1),err(1,2,2),'k.');
    legend('Experimental','Prediction','Location','northwest');

else

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
    title(sprintf(['Global Optimization for Matching Markov Model to ', ProfSel ,' Experimental Profile \n (RMSE = %0.2e, leakage = %0.2e)'],mean(errors),leakage));
    if randomModelflag
        randerr = OptimizationResults.(ProfSel).RandomResults.error  ;
        title(sprintf(['Global Optimization for Matching Markov Model to ', ProfSel ,' Experimental Profile \n (Fitted RMSE = %0.2e, Random RMSE = %0.2e, leakage = %0.2e)'],mean(errors),median(randerr),leakage));
    end
    set(gca,'fontweight','bold')
    set(gca,'TickLength',[0 0]);
end

hold off


end