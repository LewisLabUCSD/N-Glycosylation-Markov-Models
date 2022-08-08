function OptimizationResults = GenerateRandomModels(Prof,GenericNetwork,DataSet,OptimizationResults,simNum)
%% Prepare fitted variables/experimental for computing model characteristics

%%%%%%%%%%%%%%%%%%%%% Extract experimental data %%%%%%%%%%%%%%%%%%%%%
ExpData = DataSet.profiles(:,strcmp(DataSet.ProfNames,Prof));

%%%%%%%%%%%%%%%%%%%%% Initate variables %%%%%%%%%%%%%%%%%%%%%
xval = rand(simNum,size(OptimizationResults.(Prof).xval,2)).*10-5;
RandomResults = struct;
StericFlag = OptimizationResults.(Prof).OptimizationProblem.StericFlag; 
AppliedGeneidx = OptimizationResults.(Prof).OptimizationProblem.AppliedGeneidx;
stericRxns = OptimizationResults.(Prof).OptimizationProblem.stericRxns  ;
Rxn_idx = OptimizationResults.(Prof).OptimizationProblem.Rxn_idx;
LacNAcLenPenalty = rand([simNum,1])*6;

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

for a = 1:simNum
    
    [error,Predata_noRes,Predata_raw,Glys_raw,PseudoConc,Glys_conc,PseudoFlux,Rxn_flux] = ApplyTPstoGenericModels(xval(a,:),Rxn_idx,GenericNetwork,ExpData,StericFlag,AppliedGeneidx,stericRxns,LacNAcLenPenalty(a));
    
    % Initiate variables
    if a == 1
    RandomResults.error = error;RandomResults.Predata_noRes = Predata_noRes';RandomResults.Predata_raw = Predata_raw;RandomResults.Glys_raw = Glys_raw;
    RandomResults.PseudoConc = PseudoConc;RandomResults.Glys_conc = Glys_conc;RandomResults.PseudoFlux = PseudoFlux';RandomResults.Rxn_flux = Rxn_flux;
    % Update variables
    else
        RandomResults.error = [RandomResults.error;error];
        RandomResults.Predata_noRes = [RandomResults.Predata_noRes;Predata_noRes'];
        RandomResults.Predata_raw = [RandomResults.Predata_raw;Predata_raw];
        RandomResults.PseudoConc = [RandomResults.PseudoConc;PseudoConc];
        RandomResults.PseudoFlux = [RandomResults.PseudoFlux;PseudoFlux'];
    end
end

% Record experimental info
RandomResults.ExpData = ExpData';
RandomResults.AnnotatedGlycans =  DataSet.LinkageResStruct(DataSet.LinkageResStructSel(:,strcmp(Prof,DataSet.ProfNames)));
RandomResults.AnnotatedMz =  DataSet.mz(DataSet.LinkageResStructSel(:,strcmp(Prof,DataSet.ProfNames)));
RandomResults.mz_all = DataSet.mz_all; 
RandomResults.xval = xval;
RandomResults.LacNAcLenPenalty = LacNAcLenPenalty;
OptimizationResults.(Prof).RandomResults = RandomResults;
end