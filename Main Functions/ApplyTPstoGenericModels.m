function [error,Predata_noRes,Predata_raw,Glys_raw,PseudoConc,Glys_conc,PseudoFlux,Rxn_flux] = ApplyTPstoGenericModels(scales,Rxn_idx,GenericNetwork,ExpData,StericFlag,AppliedGeneidx,stericRxns)
%% Extract variables
TM = GenericNetwork.TM;
Geneidx = GenericNetwork.Geneidx;
AbsGlyIdx = GenericNetwork.AbsGlyIdx;
Glys = GenericNetwork.Glys;
pi0 = GenericNetwork.pi0;
TMidx = GenericNetwork.AllrxnList_TMidx;
[TMrow,TMcol] = ind2sub(size(TM),TMidx);
AllrxnList_steric = GenericNetwork.AllrxnList_steric;

%% Modify transition probability matrix based on the fed scaling factor

% scale rxns
for k = 1:length(Rxn_idx)
    TM(Geneidx{k}) = (10.^scales(k)).*TM(Geneidx{k});
end

% consider steric interactions if StericFlag is true
if StericFlag
    for k = length(AppliedGeneidx)-length(stericRxns)+1:length(AppliedGeneidx)
        TM(AppliedGeneidx{k}) = (1./exp(AllrxnList_steric(AppliedGeneidx{k}).*scales(k-3))).*TM(AppliedGeneidx{k});
    end
end

%% Compute stationary distribution profile from the scaled Markov model
mc = dtmc(TM);
result = redistribute(mc,25,'X0',pi0);

%% Compute & Record Model characteristics
Predata_raw = result(end,:); % Stationary distribution of all glycans, isoform-specific
Predata_noRes = cellfun(@(x) sum(Predata_raw(x)),AbsGlyIdx);% Stationary distribution of all glycans (do not distinguish isoforms)

% Stationary distribution of all glycans, isoform-specific
Glys_raw = Glys(Predata_raw~=0); % glycan labels for
Predata_raw = Predata_raw(Predata_raw~=0);

% Pseudo-concentrations of fall intermediate glycans, isoform-specific
PseudoConc = max(result); % pseudo-concentration
TM = TM./sum(TM,2);
PseudoFlux = PseudoConc'.* full(TM); % psuedo-fluxes
PseudoFlux = PseudoFlux(TMidx);
Rxn_flux = [Glys(TMrow),Glys(TMcol)];

Glys_conc = Glys(PseudoConc~=0); % glycan labels for pseudo-concentration
PseudoConc = PseudoConc(PseudoConc~=0);


%% Compute objective function
error = sum((ExpData-Predata_noRes).^2)./sqrt(length(ExpData));

end
