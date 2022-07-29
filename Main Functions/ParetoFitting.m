function errors = PSFitting(scales,TM,Geneidx,AbsGlyIdx,ExpData,pi0,mzRes,linkagePosRes)


%% Modify transition probability matrix based on the fed scaling factor

% scale rxns

for k = 1:length(Geneidx)
    TM(Geneidx{k}) = (10.^scales(k)).*TM(Geneidx{k});
end

%% Compute stationary distribution profile from the scaled Markov model
mc = dtmc(TM);
result = redistribute(mc,30,'X0',pi0);
Predata_raw = result(end,:);

Predata_noRes = cellfun(@(x) sum(Predata_raw(x)),AbsGlyIdx);
leakage = ExpData - Predata_noRes;

% impose linkage restriction
AbsGlyIdx(mzRes) = linkagePosRes;
Predata = cellfun(@(x) sum(Predata_raw(x)),AbsGlyIdx);

%% Compute objective function
error1 = (sum((ExpData-Predata_noRes).^2) + sum(leakage.^2))./sqrt(length(ExpData));
error2 = sum((Predata_noRes(Predata_noRes~=0)-Predata(Predata_noRes~=0)).^2)./sqrt(length(ExpData));
errors = [error1;error2];

end
