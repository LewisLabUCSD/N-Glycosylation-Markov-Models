function [xval_comb,allLacNAcPenalty] = PredictCombKOProfiles(KOs,BaseProfSel,ComparativeResults,OptimizationResults,GenericNetwork)
%% Source control
Enz = {'St3gal','B4galt','Mgat1','Mgat2','Mgat4','Mgat5','B3gnt'};
Rxns = {'a3SiaT','b4GalT','GnTI','GnTII','GnTIV','GnTV','iGnT'};
allEnzIdx = zeros(size(KOs));
for a = 1:length(KOs)
    allEnzIdx(a) = find(cellfun(@(x) contains(KOs{a},x),Enz),1);
end

basexval = OptimizationResults.(BaseProfSel).xval;
TPFCVec = zeros(length(KOs),size(basexval,2));
TPFCVecSpread = zeros(length(KOs),size(basexval,2));
LacNAcPenaltyFC = zeros(1,length(KOs));
Sigflag = zeros(length(KOs),size(basexval,2));
RxnTypes = GenericNetwork.RxnTypes;

%% Extract single-KO info and predict combinatorial knockouts
% Load variables
for b = 1:length(KOs)

    % load fold change data computed from Step 5
    FCInfo = ComparativeResults.([KOs{b},'_vs_',BaseProfSel]);
    RxnComp = FCInfo.RxnComp;
    [~,RxnIdx] = ismember(RxnComp,RxnTypes);
    TPFCVec(b,RxnIdx) = FCInfo.TPComp.TPFCMean;
    TPFCVecSpread(b,RxnIdx) = FCInfo.TPComp.TPFCSpread;
    LacNAcPenaltyFC(b) = FCInfo.LacNAcPenaltyComp.LacNAcPenaltyFC;
    Sigflag(b,RxnIdx) = FCInfo.OverallSig;

end
LacNAcPenaltyFC = mean(LacNAcPenaltyFC);

% Heuristically determine the new TP parameters
SampSize = 500;
FCVec_Tot = zeros(SampSize,size(basexval,2));
for k = 1:SampSize
    FCVec = zeros(1,size(basexval,2));
    for a = 1:length(RxnTypes)
        for b = 1:length(KOs)

            if Sigflag(b,a)
                EnzIdx = find(cellfun(@(x) contains(KOs{b},x) ,Enz),1);
                RxnsIdx = find(cellfun(@(x) contains(RxnTypes{a},x) ,Rxns),1);
                if EnzIdx == RxnsIdx
                    FCVec(a) = FCVec(a)+(normrnd(TPFCVec(b,a),TPFCVecSpread(b,a)));
                elseif  ~any(RxnsIdx==allEnzIdx)
                    FCVec(a) = FCVec(a)+normrnd(TPFCVec(b,a),TPFCVecSpread(b,a))./sum(Sigflag(:,a));
                end
            end

        end
    end
    FCVec_Tot(k,:) = 2.^FCVec;
end

% Apply perturbation
SampSize = 40;
RxnSize = size(basexval,2);
basevalSize = size(basexval,1);
basexval_tot = ones(SampSize,RxnSize);
allLacNAcPenalty = ones(SampSize,1);
for k = 1:SampSize
    basexval_temp = 10.^basexval(randi(basevalSize),:);
    basexval_temp = basexval_temp.*FCVec_Tot(randi(size(FCVec_Tot,1)),:);
    for a = 1:size(basexval_temp,1)
        for b = 1:size(basexval_temp,2)
            if log10(basexval_temp(a,b))>8
                basexval_temp(a,b) = 10e8;
            elseif log10(basexval_temp(a,b))<-8
                basexval_temp(a,b) = 10e-8;
            end
        end
    end
    basexval_tot(k,:) = basexval_temp;
    allLacNAcPenalty(k) = OptimizationResults.(BaseProfSel).LacNAcLenPenalty(randi(basevalSize)).*LacNAcPenaltyFC;
end
xval_comb = log10(basexval_tot);

end