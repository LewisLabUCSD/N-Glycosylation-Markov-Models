function [xval_comb,allLacNAcPenalty] = PredictCombKOProfiles(KOs,BaseProfSel,ComparativeResults,OptimizationResults,GenericNetwork)

% Root control
Enz = {'St3gal','B4galt','Mgat5','Mgat4','Mgat2','Mgat1','B3gnt'};
Rxns = {'a3SiaT','b4GalT','GnTV','GnTIV','GnTII','GnT1','iGnT'};

basexval = OptimizationResults.(BaseProfSel).xval;
basesensitivity = OptimizationResults.(BaseProfSel).SensitivityAnalysis;  
FCVec = zeros(1,size(basexval,2));
allLacNAcPenalty = [];

% Extract single-KO info and predict combinatorial knockouts
for b = 1:length(KOs)

    % load fold change data computed from Step 5
    FCInfo = ComparativeResults.([KOs{b},'_vs_',BaseProfSel]);
    RxnComp = FCInfo.RxnComp;
    FluxSig = FCInfo.PseudFluxComp.FluxFCSig;
    SensiSig = FCInfo.SensitivityComp.SensitivityMatSig; SensiSig = SensiSig(:,1) | SensiSig(:,2);
    TPFC = FCInfo.TPComp.TPFCMean;TPFCSig = FCInfo.TPComp.TPFCSig;
    allLacNAcPenalty = [allLacNAcPenalty,FCInfo.LacNAcPenaltyComp.LacNAcPenaltyFC];
    Sigflag = FluxSig & SensiSig & TPFCSig;

    for c = 1:length(RxnComp)

        RxnIdx = find(strcmp(RxnComp{c},GenericNetwork.RxnTypes),1);

        RxnTypeIdx = find(cellfun(@(x) contains(RxnComp{c},x), Rxns ),1);
        EnzIdx = find(cellfun(@(x) contains(KOs{b},x), Enz ),1);

        % if the reaction is directly related to the knockout
        if ~isempty(RxnTypeIdx) && ~isempty(EnzIdx)

            if RxnTypeIdx == EnzIdx && Sigflag(c)
                FCVec(1,RxnIdx) = FCVec(1,RxnIdx) + TPFC(c);
            elseif RxnTypeIdx ~= EnzIdx && Sigflag(c)
                newFC = TPFC(c);
                if FCVec(1,RxnIdx) == 0
                    FCVec(1,RxnIdx) = newFC;
                else
                    if sign(FCVec(1,RxnIdx)) == sign(newFC)
                        FCVec(1,RxnIdx) = sign(newFC)*max(abs([newFC,FCVec(1,RxnIdx)]));
                    else
                        FCVec(1,RxnIdx) = FCVec(1,RxnIdx)+newFC;
                    end
                end
            end


        end

    end

end

xval_comb = basexval+FCVec;
allLacNAcPenalty = max(allLacNAcPenalty);
allLacNAcPenalty = (OptimizationResults.(BaseProfSel).LacNAcLenPenalty).*allLacNAcPenalty;

end