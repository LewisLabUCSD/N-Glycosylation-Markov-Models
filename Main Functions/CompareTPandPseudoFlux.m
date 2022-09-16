function ComparativeResults =  CompareTPandPseudoFlux(Prof1, Prof2,OptimizationResults,GenericNetwork,ComparativeResults)
%% Load & initiate variables
RxnTypes = GenericNetwork.RxnTypes;
[RxnTypes,RxnIdx] = setdiff(RxnTypes,{'cg2mg','mg2tg','tg2ab'},'stable');
RxnIdx = RxnIdx';
TP1 = 10.^OptimizationResults.(Prof1).xval(:,RxnIdx);
TP2 = 10.^OptimizationResults.(Prof2).xval(:,RxnIdx);
Flux1 = OptimizationResults.(Prof1).FluxesbyComp(:,RxnIdx);
Flux2 = OptimizationResults.(Prof2).FluxesbyComp(:,RxnIdx);
Sens1 = OptimizationResults.(Prof1).SensitivityAnalysis.ErrorOverPertPct;
Sens2 = OptimizationResults.(Prof2).SensitivityAnalysis.ErrorOverPertPct;
Sens1Rxns = OptimizationResults.(Prof1).SensitivityAnalysis.RxnNames;
Sens2Rxns = OptimizationResults.(Prof2).SensitivityAnalysis.RxnNames;
perturbationPct = OptimizationResults.(Prof2).SensitivityAnalysis.perturbationPct;

TPFCMean = zeros(size(RxnTypes));
TPFCSig = zeros(size(RxnTypes));
FluxFCMean = zeros(size(RxnTypes));
FluxFCSig = zeros(size(RxnTypes));
SensitivityMat = zeros(length(RxnTypes),2);
SensitivityMatSig = SensitivityMat;

% Set significance level
pval = 0.05;
sensSig = 0.1;

%% Compare model TPs

for a = 1:length(RxnTypes)


    % compute sensitivity significance
    indexSel1 = find(strcmp(RxnTypes{a},Sens1Rxns),1);
    indexSel2 = find(strcmp(RxnTypes{a},Sens2Rxns),1);
    Sens1_temp = abs(Sens1(indexSel1,:));Sens2_temp = abs(Sens2(indexSel2,:));
    SensitivityMat(a,:) = [max(Sens1_temp),max(Sens2_temp)];
    SensitivityMatSig(a,:) = SensitivityMat(a,:)>sensSig;

    % Compute TP fold change
    TP1Min = min(perturbationPct(Sens1_temp<0.01));TP1Max = max(perturbationPct(Sens1_temp<0.01));
    if isempty(TP1Min);TP1Min = 0;end
    if isempty(TP1Max);TP1Max = 0;end
    if any([TP1Min~=0,TP1Max~=0])
        TP1PerturbRange = linspace(min(perturbationPct(Sens1_temp<0.01)),max(perturbationPct(Sens1_temp<0.01)),50);
        TP1_temp = zeros(length(TP1(:,a))*50,1);
        for k = 1:50
            TP1_temp((k-1)*length(TP1(:,a))+1:k*length(TP1(:,a))) = log10(TP1(:,a)).*(1+TP1PerturbRange(k));
        end
        % TP1_temp = 10.^TP1_temp;
    else
        TP1_temp = log10(TP1(:,a));
    end

    TP2Min = min(perturbationPct(Sens2_temp<0.01));TP2Max = max(perturbationPct(Sens2_temp<0.01));
    if isempty(TP2Min);TP2Min = 0;end
    if isempty(TP2Max);TP2Max = 0;end
    if any([TP2Min~=0,TP2Max~=0])
        TP2PerturbRange = linspace(min(perturbationPct(Sens2_temp<0.01)),max(perturbationPct(Sens2_temp<0.01)),50);
        TP2_temp = zeros(length(TP2(:,a))*50,1);
        for k = 1:50
            TP2_temp((k-1)*length(TP2(:,a))+1:k*length(TP2(:,a))) = log10(TP2(:,a)).*(1+TP2PerturbRange(k));
        end
        % TP2_temp = 10.^TP2_temp;
    else
        TP2_temp = log10(TP2(:,a));
    end
    
    % TP1_temp = TP1(:,a);TP2_temp = TP2(:,a);
    TPFCDist = TP1_temp-TP2_temp';TPFCDist = TPFCDist(:);
    TPCombVec = [TP1_temp;TP2_temp];TPCombVecLen = length(TPCombVec);
    TPFCPermDistq = zeros(1000,1);
    for b = 1:1000 % 1000 permutations
        TPFCPermDistq(b) = TPCombVec(randi(TPCombVecLen))-TPCombVec(randi(TPCombVecLen));
    end
    TPFCMean(a) = median(TPFCDist);
    TPFCSig(a) = ranksum(TPFCPermDistq,TPFCDist)<pval;

    % Comptue Flux fold change
    FluxFCMean(a) = log10(median(Flux1(:,a))/median(Flux2(:,a)));
    FluxFCDist = Flux1(:,a)./Flux2(:,a)';FluxFCDist = FluxFCDist(:);
    FluxCombVec = [Flux1(:,a);Flux2(:,a)];FluxCombVecLen = length(FluxCombVec);
    FluxFCPermDistq = zeros(1000,1);
    for b = 1:1000 % 1000 permutations
        FluxFCPermDistq(b) = FluxCombVec(randi(FluxCombVecLen))./FluxCombVec(randi(FluxCombVecLen));
    end
    FluxFCMean(a) = log10(median(FluxFCDist));
    FluxFCSig(a) = ranksum(FluxFCPermDistq,FluxFCDist)<pval;

end

%% Compare LacNAc elongation factor
if ~isempty(OptimizationResults.(Prof1).LacNAcLenPenalty)
    Pen1 = OptimizationResults.(Prof1).LacNAcLenPenalty;
    Pen2 = OptimizationResults.(Prof2).LacNAcLenPenalty;

    PenSens1 = OptimizationResults.(Prof1).SensitivityAnalysis.ErrorOverPertPct(strcmp(OptimizationResults.(Prof1).SensitivityAnalysis.RxnNames,'LacNAc Penalty'),:);
    PenSens2 = OptimizationResults.(Prof2).SensitivityAnalysis.ErrorOverPertPct(strcmp(OptimizationResults.(Prof2).SensitivityAnalysis.RxnNames,'LacNAc Penalty'),:);

    flag1 = ranksum(Pen1,Pen2)<pval;

    if flag1 && max(abs(PenSens1)) < sensSig && max(abs(PenSens2)) > sensSig
        LacNAcPenaltyFC = median(PenSens1);
    elseif flag1 && max(abs(PenSens1)) > sensSig && max(abs(PenSens2)) < sensSig
        LacNAcPenaltyFC = median(PenSens2);
    elseif flag1 && max(abs(PenSens1)) < sensSig && max(abs(PenSens2)) < sensSig
        LacNAcPenaltyFC = max([PenSens1,PenSens2]);
    else
        LacNAcPenaltyFC = 1;
    end
end

%% Plot comparison results

Prof1 = strrep(Prof1,'_','/');
Prof2 = strrep(Prof2,'_','/');

% TP
figure;
hold on
h = subplot(1,3,1);
hold on
imagesc(TPFCMean);
colormap redbluecmap
climit = max(abs(TPFCMean));
caxis([-climit climit]);
colorbar;
yticks(1:length(RxnTypes));
yticklabels(strrep(RxnTypes,'_',' '));
xticks([]);
for a = 1:length(TPFCSig)
    if TPFCSig(a) && (TPFCMean(a)>log10(1.1) || TPFCMean(a)<log10(0.9))
        plot(1,a,'o','MarkerFaceColor','#f7c602','MarkerEdgeColor','#f7c602');
    end
end
title(['log10(TP_{',Prof1,'}-log10(TP_{',Prof2,'})']);
ylabel('Reaction Types','FontWeight','bold');
hold off

% Flux
subplot(1,3,2)
hold on
imagesc(FluxFCMean);
colormap redbluecmap
climit = max(abs(FluxFCMean));
caxis([-climit climit]);
colorbar;
yticks([]);
xticks([]);
for a = 1:length(FluxFCSig)
    if FluxFCSig(a) && (FluxFCMean(a)>log10(1.1) || FluxFCMean(a)<log10(0.9))
        plot(1,a,'o','MarkerFaceColor','#f7c602','MarkerEdgeColor','#f7c602');
    end
end
title(['log10(Flux_{',Prof1,'}/Flux_{',Prof2,'})']);
hold off

% Sensitivity
ax = subplot(1,3,3);
hold on
imagesc(SensitivityMat);
colormap(ax, 'autumn');
climit = abs(SensitivityMat(:));
climit = max(climit(~isoutlier(climit)));
caxis([0 climit]);
colorbar;
yticks([]);
for a = 1:size(SensitivityMat,1)
    for b = 1:size(SensitivityMat,2)
        if SensitivityMatSig(a,b)
            plot(b,a,'o','MarkerFaceColor','#f7c602','MarkerEdgeColor','#f7c602');
        end
    end
end
xticks(1:2);
xticklabels({Prof1,Prof2});
xlabel('Prof Name','FontWeight','bold');
title({'Sensitivity for each reaction','(max abs(RMSE%/Perturbation%))'});
hold off
hold off

%% Highlight significant reaction types
ticklabels = get(h,'YTickLabel');
OverallSig = false(length(ticklabels),1);
ticklabels_new = cell(size(ticklabels));
for k = 1:length(ticklabels)
    if (TPFCSig(k) && (TPFCMean(k)>log10(1.1) || TPFCMean(k)<log10(0.9))) ...
                 && (FluxFCSig(k) && (FluxFCMean(k)>log10(1.1) || FluxFCMean(k)<log10(0.9)))...
              && any(SensitivityMatSig(k,:))

        ticklabels_new{k} = ['\color{red} ' ticklabels{k}];
        OverallSig(k) = true;
    else
        ticklabels_new{k} = ticklabels{k};
    end
end
set(h, 'YTickLabel', ticklabels_new);

%% Record results
CompName = strrep([Prof1 '_vs_' Prof2],'/','_');
ComparativeResults.(CompName).TPComp.TPFCMean = TPFCMean;
ComparativeResults.(CompName).TPComp.TPFCSig = TPFCSig;
ComparativeResults.(CompName).PseudFluxComp.FluxFCMean = FluxFCMean;
ComparativeResults.(CompName).PseudFluxComp.FluxFCSig = FluxFCSig;
ComparativeResults.(CompName).SensitivityComp.SensitivityMat = SensitivityMat;
ComparativeResults.(CompName).SensitivityComp.SensitivityMatSig = SensitivityMatSig;
ComparativeResults.(CompName).RxnComp = RxnTypes;
ComparativeResults.(CompName).LacNAcPenaltyComp.LacNAcPenaltyFC = LacNAcPenaltyFC;
ComparativeResults.(CompName).OverallSig = OverallSig;
end