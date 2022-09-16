function [xval_comb,allLacNAcPenalty] = PredictCombKOProfiles(KOs,BaseProfSel,ComparativeResults,OptimizationResults,GenericNetwork)

% Root control
Enz = {'St3gal','B4galt','Mgat5','Mgat4','Mgat2','Mgat1','B3gnt'};
Rxns = {'a3SiaT','b4GalT','GnTV','GnTIV','GnTII','GnTI','iGnT'};

basexval = OptimizationResults.(BaseProfSel).xval;
basesensitivity = OptimizationResults.(BaseProfSel).SensitivityAnalysis;  
FCVec = zeros(1,size(basexval,2));
LacNAcPenaltyFC = [];
Sigflag_all = [];

% Extract single-KO info and predict combinatorial knockouts
for b = 1:length(KOs)
    Sigflag_all = [Sigflag_all,ComparativeResults.([KOs{b},'_vs_',BaseProfSel]).OverallSig];
end

for b = 1:length(KOs)

    % load fold change data computed from Step 5
    FCInfo = ComparativeResults.([KOs{b},'_vs_',BaseProfSel]);
    RxnComp = FCInfo.RxnComp;
    TPFC = FCInfo.TPComp.TPFCMean;TPFCSig = FCInfo.TPComp.TPFCSig;
    LacNAcPenaltyFC = [LacNAcPenaltyFC,FCInfo.LacNAcPenaltyComp.LacNAcPenaltyFC];
    Sigflag = FCInfo.OverallSig;

    for c = 1:length(RxnComp)

        RxnIdx = find(strcmp(RxnComp{c},GenericNetwork.RxnTypes),1);
        RxnTypeIdx = find(cellfun(@(x) contains(RxnComp{c},x), Rxns ),1);
        EnzIdx = find(cellfun(@(x) contains(KOs{b},x), Enz ),1);
        
        if ~isempty(RxnTypeIdx) && ~isempty(EnzIdx)
            
            % if the reaction is directly related to the knockout
            if RxnTypeIdx == EnzIdx && Sigflag(c)
                FCVec(1,RxnIdx) = FCVec(1,RxnIdx) + TPFC(c);
            % if the reaction is indirectly related to the knockout
            elseif RxnTypeIdx ~= EnzIdx && Sigflag(c)                
                FCVec(1,RxnIdx) = FCVec(1,RxnIdx)+ TPFC(c)/sum(Sigflag_all(c,:));
            end
        end
    end
end

% Apply perturbation
basexval = 10.^basexval;
basexval = basexval.*10.^FCVec;
for a = 1:size(basexval,1)
    for b = 1:size(basexval,2)
        if log10(basexval(a,b))>8
           basexval(a,b) = 10e8;
        elseif log10(basexval(a,b))<-8
            basexval(a,b) = 10e-8;
        end
    end
end

xval_comb = log10(basexval);
LacNAcPenaltyFC = min(LacNAcPenaltyFC);
allLacNAcPenalty = (OptimizationResults.(BaseProfSel).LacNAcLenPenalty).*LacNAcPenaltyFC;

end