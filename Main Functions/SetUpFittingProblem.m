function OptimizationProblem = SetUpFittingProblem(num,Prof,GenericNetwork,DataSet,StericFlag,UseWTStericFlag,WTSteric)
%% Prepare fitting variables

%%%%%%%%%%%%%%%%%%%%% Load data & normalize experimental data (sum of signal)%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% intensities in each % profile is equal to 1) %%%%%%%%%%%%%%%%%%%%%
ExpData = DataSet.profiles(:,strcmp(DataSet.ProfNames,Prof));
ExpData = ExpData./sum(ExpData);
Geneidx = GenericNetwork.Geneidx;
TM = GenericNetwork.TM;
AbsGlyIdx = GenericNetwork.AbsGlyIdx;
pi0 = GenericNetwork.pi0;
AllrxnList_steric = GenericNetwork.AllrxnList_steric;
stericRxns = GenericNetwork.stericRxns;

%%%%%%%%%%%%%%%%%%%%% Initiate temporary variables for storing fitting results %%%%%%%%%%%%%%%%%%%%%
% steric interactions
for a = 1:length(stericRxns)
    Geneidx = [Geneidx;{find(contains(GenericNetwork.AllrxnList_RxnTypes,stericRxns{a}))}];
end

% Modify Geneidx to include steric parameters
if StericFlag && ~UseWTStericFlag
    xval = zeros(num,length(Geneidx));
else
    xval = zeros(num,length(Geneidx)-length(stericRxns));
end
fval = zeros(num,1);

%% Prepare fitting constraints and parameters

% Define optimization constraints (refer to the help document for MATLAB
% function "optimproblem"). The scaling factor for each reaction type
% (RxyTypes) is in the range of 10^[-4,4].
Rxn_idx = find(~ismember(GenericNetwork.RxnTypes,{'cg2mg','mg2tg','tg2ab'}));
optimproblem.x0 = 10*rand(1,size(xval,2)-3)-5;
optimproblem.lb = ones(1,size(xval,2)-3).*-5;
if StericFlag && ~UseWTStericFlag
    optimproblem.lb(end-length(stericRxns)+1:end) = 0;
end
optimproblem.ub = ones(1,size(xval,2)-3).*5;
optimproblem.Aineq = [];
optimproblem.bineq = [];
optimproblem.Aeq = [];
optimproblem.beq = [];
optimproblem.nonlcon = [];
optimproblem.intcon = [];
optimproblem.rngstate = [];

% Define optimization parameters (refer to the help document for MATLAB
% function "patternsearch").
optimproblem.solver = 'patternsearch'; % choose pattern search algorithm
optimproblem.options = optimoptions('patternsearch',...
    'MaxTime',7200,...
    'Display','iter',...
    'UseParallel',true,...
    'PollOrderAlgorithm','random',...
    'MeshRotate','on',...
    'MeshContractionFactor',0.75,...
    'MeshExpansionFactor',2,...
    'AccelerateMesh',true,...
    'ScaleMesh',false,...
    'UseCompletePoll',true,...
    'UseCompleteSearch',true,...
    'InitialMeshSize',1e8,...
    'Cache','off',...
    'FunctionTolerance',1e-9,...
    'MaxIterations',8000,...
    'MaxFunctionEvaluations',8000*length(optimproblem.x0),...
    'SearchFcn','GSSPositiveBasis2N',... % GSSPositiveBasisNp1 MADSPositiveBasis2N MADSPositiveBasisNp1
    'UseParallel',true,...
    'PlotFcn',{'psplotbestf','psplotbestx'});

%% Modify AbsGlyIdx based on glycan annotations
%[mzRes, linkagePosRes] = PrepGlycanAnnotationConstraints(Prof, DataSet, GenericNetwork);

% Input:
% 1. Prof: the string of the selected profile name for fitting (must also be an element in DataSet.ProfNames
% 2. DataSet: struct variable contructed from Step 1.
% 3. GenericNetwork: struct variable constructed from Step 2.

% Output:
% 1. mzRes: double-type indices of m/z values at which glycan structures were annotated
% 2. linkagePosRes: cell of strings representing the linear codes of actual glycan annotations
% corresponding to the m/z values (mzRes)

[mzRes, linkagePosRes] = PrepGlycanAnnotationConstraints(Prof, DataSet, GenericNetwork);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up objective function
if ~UseWTStericFlag
    optimproblem.objective  = @(x) PSFitting(x,TM,Geneidx,Rxn_idx,AbsGlyIdx,ExpData,pi0,mzRes,linkagePosRes,StericFlag,AllrxnList_steric);
else
    optimproblem.objective  = @(x) PSFitting(x,TM,Geneidx,Rxn_idx,AbsGlyIdx,ExpData,pi0,mzRes,linkagePosRes,StericFlag,AllrxnList_steric,WTSteric);
end

%% Record the optimization problem setting
OptimizationProblem.optimproblem = optimproblem;
OptimizationProblem.mzRes = mzRes;
OptimizationProblem.linkagePosRes = linkagePosRes;
OptimizationProblem.NumOfModels = num;
OptimizationProblem.ExpData = ExpData;
OptimizationProblem.Rxn_idx = Rxn_idx;
OptimizationProblem.stericRxns = stericRxns;
OptimizationProblem.AppliedGeneidx = Geneidx;
OptimizationProblem.StericFlag = StericFlag;
OptimizationProblem.UseWTStericFlag = UseWTStericFlag;
OptimizationProblem.WTSteric = WTSteric;
OptimizationProblem.xval = xval;
OptimizationProblem.fval = fval;

end