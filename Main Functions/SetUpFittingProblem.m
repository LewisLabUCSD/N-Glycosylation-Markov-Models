function OptimizationProblem = SetUpFittingProblem(num,Prof,GenericNetwork,DataSet,StericFlag,UseWTStericFlag,WTSteric,Method)
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
AllrxnList_LacNAcLen = GenericNetwork.AllrxnList_LacNAcLen;
AllrxnList_LacNAcLen_idx = GenericNetwork.AllrxnList_LacNAcLen_idx;
LeakageGlyIdx = GenericNetwork.LeakageGlyIdx;

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
optimproblem.x0 = 16*rand(1,size(xval,2)-3+1)-8;
% optimproblem.x0 = ones(1,size(xval,2)-3+1).*8;
optimproblem.lb = ones(1,size(xval,2)-3+1).*-8;
if StericFlag && ~UseWTStericFlag
    optimproblem.lb(end-length(stericRxns)+1:end-1) = 0;
end
optimproblem.lb(end) = 0;
optimproblem.ub = ones(1,size(xval,2)-3+1).*8;
optimproblem.Aineq = [];
optimproblem.bineq = [];
optimproblem.Aeq = [];
optimproblem.beq = [];
optimproblem.nonlcon = [];
optimproblem.intcon = [];
optimproblem.rngstate = [];

% Define optimization parameters (refer to the help document for MATLAB
% function 'patternsearch' and 'simulannealbnd').

if strcmp(Method,'PatternSearch')
    
    optimproblem.solver = 'patternsearch'; % choose pattern search algorithm
    optimproblem.options = optimoptions('patternsearch',...
        'MaxTime',7200,...
        'Display','iter',...
        'PollingOrder','Random',...
        'MeshRotate','on',...
        'MeshContraction',0.75,...
        'MeshExpansionFactor',3,...
        'AccelerateMesh',false,...
        'ScaleMesh',true,...
        'UseCompletePoll',true,...
        'UseCompleteSearch',true,...
        'InitialMeshSize',1e8,...
        'Cache','off',...
        'FunctionTolerance',1e-7,...
        'MaxIterations',7000,...
        'MaxFunctionEvaluations',20000*length(optimproblem.x0),...
        'SearchFcn','MADSPositiveBasis2N',... % GSSPositiveBasisNp1 MADSPositiveBasis2N MADSPositiveBasisNp1 searchlhs searchga searchneldermead
        'PollMethod','GSSPositiveBasis2N',...
        'UseParallel',true,...
        'PlotFcn',{'psplotbestf','psplotbestx'});

elseif strcmp(Method,'PSPSHybrid')
    
    optimproblem.nvars = size(xval,2)-3;
    hybridoptions = optimoptions('patternsearch',...
        'MaxTime',4000,...
        'Display','iter',...
        'PollingOrder','Random',...
        'MeshRotate','on',...
        'MeshContraction',0.75,...
        'MeshExpansionFactor',2,...
        'AccelerateMesh',true,...
        'ScaleMesh',true,...
        'UseCompletePoll',true,...
        'UseCompleteSearch',true,...
        'InitialMeshSize',1e8,...
        'Cache','off',...
        'FunctionTolerance',1e-7,...
        'MaxIterations',5000,...
        'MaxFunctionEvaluations',8000*length(optimproblem.x0),...
        'SearchFcn','GSSPositiveBasis2N',... % GSSPositiveBasisNp1 MADSPositiveBasis2N MADSPositiveBasisNp1 searchlhs searchga searchneldermead
        'PollMethod','MADSPositiveBasis2N',...
        'UseParallel',true,...
        'PlotFcn',{'psplotbestf','psplotbestx'});

    optimproblem.solver = 'particleswarm'; % choose pattern search algorithm
    optimproblem.options = optimoptions('particleswarm',...
        'MaxTime',3200,...
        'Display','iter',...
        'FunctionTolerance',1e-7,...
        'MaxIterations',5000,...
        'HybridFcn',{@patternsearch, hybridoptions},...
        'MaxStallTime',600,...
        'PlotFcn',{'pswplotbestf'},...
        'SwarmSize',100,...
        'UseParallel',true);
end

%% Read glycan annotations
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
optimproblem.objective  = @(x) PSFitting(x,TM,Geneidx,Rxn_idx,AbsGlyIdx,LeakageGlyIdx,ExpData,pi0,mzRes,linkagePosRes,StericFlag,AllrxnList_steric,WTSteric,AllrxnList_LacNAcLen_idx,AllrxnList_LacNAcLen);

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