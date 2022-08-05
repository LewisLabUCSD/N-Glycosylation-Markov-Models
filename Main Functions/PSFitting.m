function error = PSFitting(scales,TM,Geneidx,Rxn_idx,AbsGlyIdx,ExpData,pi0,mzRes,linkagePosRes,StericFlag,AllrxnList_steric,WTSteric,AllrxnList_LacNAcLen_idx,AllrxnList_LacNAcLen)
%% Modify transition probability matrix based on the fed scaling factor

% scale rxns
for k = 1:length(Rxn_idx)
    TM(Geneidx{Rxn_idx(k)}) = (10.^scales(k)).*TM(Geneidx{Rxn_idx(k)});
end

% Consider LacNAc length (the last variable)
TM(AllrxnList_LacNAcLen_idx)  = TM(AllrxnList_LacNAcLen_idx)./exp(scales(end).*AllrxnList_LacNAcLen);

%% Consider steric interactions 
if StericFlag && isempty(WTSteric)
    for k = length(Rxn_idx)+4:length(Geneidx)
        TM(Geneidx{k}) = (1./exp(AllrxnList_steric(Geneidx{k}).*scales(k-3))).*TM(Geneidx{k});
    end
end

% apply existing steric factors from fitted WT models
if StericFlag && ~isempty(WTSteric)
    count = 1;
    for k = length(Rxn_idx)+4:length(Geneidx)
        TM(Geneidx{k}) = (1./exp(AllrxnList_steric(Geneidx{k}).*WTSteric(count))).*TM(Geneidx{k});
        count = 1+count;
    end
end

%% Compute stationary distribution profile from the scaled Markov model
mc = dtmc(TM);
result = redistribute(mc,25,'X0',pi0);
Predata_raw = result(end,:);

Predata_noRes = cellfun(@(x) sum(Predata_raw(x)),AbsGlyIdx);

% impose linkage restriction
if ~isempty(mzRes)
    AbsGlyIdx(mzRes) = linkagePosRes;
    Predata = cellfun(@(x) sum(Predata_raw(x)),AbsGlyIdx);
end

%% Process leakage signals
leakage = 1-sum(Predata_noRes);

%% Compute objective function
if ~isempty(mzRes)
    error = (sum((ExpData-Predata_noRes).^2) + sum((ExpData(mzRes)-Predata(mzRes)).^2)+ leakage.^2)./sqrt(length(ExpData));
else
    error = (sum((ExpData-Predata_noRes).^2) + leakage.^2)./sqrt(length(ExpData));
end

end
