function OptimizationResults = ComputeFittedModels(Prof,GenericNetwork,DataSet,OptimizationResults)
%% Prepare fitted variables/experimental for computing model characteristics

%%%%%%%%%%%%%%%%%%%%% Extract experimental data %%%%%%%%%%%%%%%%%%%%%
ExpData = DataSet.profiles(:,strcmp(DataSet.ProfNames,Prof));

%%%%%%%%%%%%%%%%%%%%% Eliminate models with outlier errors %%%%%%%%%%%%%%%%%%%%%
xval = OptimizationResults.(Prof).xval;
LacNAcLenPenalty = OptimizationResults.(Prof).LacNAcLenPenalty;
StericFlag = OptimizationResults.(Prof).OptimizationProblem.StericFlag; 
AppliedGeneidx = OptimizationResults.(Prof).OptimizationProblem.AppliedGeneidx;
stericRxns = OptimizationResults.(Prof).OptimizationProblem.stericRxns  ;
Rxn_idx = OptimizationResults.(Prof).OptimizationProblem.Rxn_idx;
AllrxnList_LacNAcLen = GenericNetwork.AllrxnList_LacNAcLen;
AllrxnList_LacNAcLen_idx = GenericNetwork.AllrxnList_LacNAcLen_idx;

%% Apply each set of fitted transition probabilities to compute model characteristics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [error,Predata_noRes,Predata_raw,Glys_raw,PseudoConc,Glys_conc,PseudoFlux,Rxn_flux] = ApplyTPstoGenericModels(xval,GenericNetwork,ExpData);

% Input:
% 1. xval: double-type vector containing the set of transition
% probabilities for a model. Extract from OptimizationResults (Step 3). 
% 2. GenericNetwork: struct  constructed from Step 2 
% 3. ExpData: a double-type column vector extracted from DataSet (Step 1),
% representing the experimental glycoprofile used for fitting

% Output:
% 1. error: the error computed by the objective function
% 2. Predata_noRes: predicted signal corresponding to m/z values in mz_all
% 3. Predata_raw: predicted signal corresponding to the absorbed glycan isoforms in Glys
% 4. Glys_raw: the absorbed glycan isoforms in Glys
% 5. PseudoConc: pseudo-concentrations of the intermediate glycans
% 6. Glys_conc: the glycan isoforms corresponding to PseudoConc entries
% 7. PseudoFlux: pseudo-fluxes of network reactions
% 8. Rxn_flux: the edge list (reactant, product) representing the network
% reactions corresponding to the fluxes in PseudoFlux

TM_avg = 0;

for a = 1:size(xval,1)
    
    [error,Predata_noRes,Predata_raw,Glys_raw,PseudoConc,Glys_conc,PseudoFlux,Rxn_flux,TM] = ApplyTPstoGenericModels(xval(a,:),Rxn_idx,GenericNetwork,ExpData,StericFlag,AppliedGeneidx,stericRxns,AllrxnList_LacNAcLen,AllrxnList_LacNAcLen_idx,LacNAcLenPenalty(a));
    TM_avg = TM_avg+TM;

    % Initiate variables
    if a == 1
    OptimizationResults.(Prof).error = error;OptimizationResults.(Prof).Predata_noRes = Predata_noRes';OptimizationResults.(Prof).Predata_raw = Predata_raw;OptimizationResults.(Prof).Glys_raw = Glys_raw;
    OptimizationResults.(Prof).PseudoConc = PseudoConc;OptimizationResults.(Prof).Glys_conc = Glys_conc;OptimizationResults.(Prof).PseudoFlux = PseudoFlux';OptimizationResults.(Prof).Rxn_flux = Rxn_flux;
    % Update variables
    else
        OptimizationResults.(Prof).error = [OptimizationResults.(Prof).error;error];
        OptimizationResults.(Prof).Predata_noRes = [OptimizationResults.(Prof).Predata_noRes;Predata_noRes'];
        OptimizationResults.(Prof).Predata_raw = [OptimizationResults.(Prof).Predata_raw;Predata_raw];
        OptimizationResults.(Prof).PseudoConc = [OptimizationResults.(Prof).PseudoConc;PseudoConc];
        OptimizationResults.(Prof).PseudoFlux = [OptimizationResults.(Prof).PseudoFlux;PseudoFlux'];
    end
end

TM_avg = TM_avg./size(xval,1);

% Record experimental info
OptimizationResults.(Prof).ExpData = ExpData';
OptimizationResults.(Prof).AnnotatedGlycans =  DataSet.LinkageResStruct(DataSet.LinkageResStructSel(:,strcmp(Prof,DataSet.ProfNames)));
OptimizationResults.(Prof).AnnotatedMz =  DataSet.mz(DataSet.LinkageResStructSel(:,strcmp(Prof,DataSet.ProfNames)));
OptimizationResults.(Prof).mz_all = DataSet.mz_all; 
OptimizationResults.(Prof).TM_avg = TM_avg;

end