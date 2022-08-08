function SensitivityAnalysis = ConductSensitivityAnalysis(Prof,GenericNetwork,DataSet,OptimizationResults,UseNumofSamples)
%% Prepare fitted variables/experimental for computing model characteristics

%%%%%%%%%%%%%%%%%%%%% Extract experimental data %%%%%%%%%%%%%%%%%%%%%
ExpData = DataSet.profiles(:,strcmp(DataSet.ProfNames,Prof));

%%%%%%%%%%%%%%%%%%%%% Eliminate models with outlier errors %%%%%%%%%%%%%%%%%%%%%
flag = ~isoutlier(OptimizationResults.(Prof).fval);
xval = OptimizationResults.(Prof).xval(flag,:);
StericFlag = OptimizationResults.(Prof).OptimizationProblem.StericFlag;
AppliedGeneidx = OptimizationResults.(Prof).OptimizationProblem.AppliedGeneidx;
stericRxns = OptimizationResults.(Prof).OptimizationProblem.stericRxns  ;
Rxn_idx = OptimizationResults.(Prof).OptimizationProblem.Rxn_idx;
LacNAcLenPenalty = OptimizationResults.(Prof).LacNAcLenPenalty;
if nargin == 5
    randSel = randi(size(xval,1),1,UseNumofSamples);
    UseNumofSamples = min([UseNumofSamples size(xval,1)]);
    xval = xval(randSel,:);
    LacNAcLenPenalty = LacNAcLenPenalty(randSel,:);
end

%% Defined perturbation percentage
% only perturb actual reactions
RxnNames = GenericNetwork.RxnTypes(Rxn_idx);
RxnNames = [RxnNames;{'LacNAc Penalty'}];
ValueRange = [-0.5 -0.1 -0.01 0 0.01 0.1 0.5];
errorMat = zeros(length(RxnNames),length(ValueRange));

f1 = waitbar(0,'Computing perturbed profile by perturbing:');
for k1 = 1:length(RxnNames)
    waitbar(k1/length(RxnNames),f1,['Computing perturbed profile by perturbing: ',RxnNames{k1}]);
    for k2 = 1:length(ValueRange)
        
        xval_temp = [xval,LacNAcLenPenalty];
        if strcmp(RxnNames{k1},'LacNAc Penalty')
            xval_temp(:,end) = xval_temp(:,end).*(1+ValueRange(k2));
        else
            xval_temp(:,Rxn_idx(k1)) = xval_temp(:,Rxn_idx(k1)).*(1+ValueRange(k2));
        end
        errorVec = zeros(size(xval_temp,1),1);
        
        parfor a = 1:size(xval_temp,1)           
            [errorVec(a)] = ApplyTPstoGenericModels(xval_temp(a,1:end-1),Rxn_idx,GenericNetwork,ExpData,StericFlag,AppliedGeneidx,stericRxns,xval_temp(a,end));
        end
        
        errorMat(k1,k2) = mean(errorVec);
        
    end
end
delete(f1);

% Normalize errorMat
errorMat = errorMat./(errorMat(:,4));
errorMat = abs(errorMat-errorMat(:,4));
errorMat = errorMat./ValueRange;
middleidx = ceil(length(ValueRange)/2);
errorMat = errorMat(:,[1:middleidx-1,middleidx+1:end]);
% sort errorMat
[~,idx] = sort(max(abs(errorMat),[],2),'descend');
errorMat = errorMat(idx,:);
RxnNames = RxnNames(idx);

%% Plot perturbation results
figure;
hold on
% plot negative
for a = 1:middleidx-1
barh(errorMat(:,a));
end
% plot positive
for a = size(errorMat,2):-1:middleidx
    barh(errorMat(:,a));
end
legend(cellstr(strcat(num2str(ValueRange([1:middleidx-1,end:-1:middleidx+1])'.*100),'%')));

% label figure
xlabel('RMSE Change%/Perturbation% (+/-)','FontWeight','bold');
ylabel('Rxn Types','FontWeight','bold');
yticks(1:length(RxnNames));
yticklabels(strrep(RxnNames,'_',' '));
title(['Sensitivity analysis for ',Prof]);
hold off

% Record data
SensitivityAnalysis.ErrorOverPertPct = errorMat;
SensitivityAnalysis.RxnNames = RxnNames;
SensitivityAnalysis.perturbationPct = ValueRange([1:middleidx-1,middleidx:end]);
end