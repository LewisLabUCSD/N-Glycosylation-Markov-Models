%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data');
load Data.mat
load GenericNetwork.mat

%% Step 3a. Define fitting parameters for Pattern Search Algorithm (stochastic global optimization)
% Please refer to the function descriptions of patternsearch and optimoptions 
% for detailed explanations. 

% Specify the names of glycoprofiles to be fitted
ProfSel  = {'WT'}; % a cell of strings of profile names selected to be fitted.
% The name must be present in GenericNetwork.ProfNames.

num = 2; % Number of models fitted for each profile
method = 'StepProgression'; % the method used to compute the stationary state of Markov models
[~,transport_idx] = ismember({'cg2mg','mg2tg','tg2ab'},GenericNetwork.RxnTypes); % indices of transports in RxnTypes
optimproblem.x0 = 10*rand(1,length(GenericNetwork.Geneidx))-5;optimproblem.x0(transport_idx) = 0; % initial point a random vector with each element in the range of [-5,5]; set transport reactions to 0
        optimproblem.lb = ones(1,length(GenericNetwork.Geneidx)).*-5;optimproblem.lb(transport_idx) = 0; % lower bound for each parameter; set transport reactions to 0
        optimproblem.ub = ones(1,length(GenericNetwork.Geneidx)).*5;optimproblem.ub(transport_idx) = 0; % upper bound for each parameter; set transport reactions to 0
        optimproblem.solver = 'patternsearch'; % choose pattern search algorithm
        
        % No other constraints
        optimproblem.Aineq = [];
        optimproblem.bineq = [];
        optimproblem.Aeq = [];
        optimproblem.beq = [];
        optimproblem.nonlcon = [];
        optimproblem.intcon = [];
        optimproblem.rngstate = [];
        
%% Step 3b. Fit Markov models to glycoprofiles by Pattern Search Algorithm (stochastic global optimization)
% Each selected glycoprofiles is fitted sequentially

for a = 1:length(ProfSel)
    
    %%%%%%%%%%%%%%%%%%%%% Normalize experimental data (sum of signal %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% intensities in each % profile is equal to 100) %%%%%%%%%%%%%%%%%%%%%
    New_ExpData = profiles(:,strcmp(ProfNames,ProfSel{a}));
    New_ExpData = New_ExpData./sum(New_ExpData).*100;
    
    %%%%%%%%%%%%%%%%%%%%% Initiate temporary variables for fitting results %%%%%%%%%%%%%%%%%%%%%
    xval = zeros(num,length(GenericNetwork.Geneidx)); % fitted transition probabilities for each "reactions" in RxnTypes
    fval = zeros(num,1); % fitting errors
    
    %%%%%%%%%%%%%%%%%%%%% Modify AbsGlyIdx based on glycan annotations %%%%%%%%%%%%%%%%%%%%%
    
    
    Annotation_idx = find(LinkageResStructSel(:,strcmp(ProfNames,ProfSel{a})));
    mzRes = mz();
    linkageRes = Glycans(Profiles(:,strcmp(KO{a},KOs))~=0);
    if strcmp(method,'StepProgression')
        linkagePosRes = cellfun(@(x) find(strcmp(x,Glys)),linkageRes,'UniformOutput',0);
    else
        linkagePosRes = cellfun(@(x) find(strcmp(x,AbsGlys),1),linkageRes,'UniformOutput',0);
    end
    
    for k = 1:num
        
        % Set up problem
        f =@(x) PSFitting(x,TM,Geneidx,AbsGlyIdx,New_ExpData,NA,NT,pi0,pi0_T,...
            AllrxnList_iGnT,AllrxnList_iGnTLen,method,mzRes,linkagePosRes);
        
        optimproblem.objective  = f;
       
        
        % Optimization options
        
                optimproblem.options = optimoptions('patternsearch',...
                    'MaxTime',1800,...
                    'Display','iter',...
                    'UseParallel',true,...
                    'PollOrderAlgorithm','random',...
                    'InitialMeshSize',1e8,...
                    'Cache','off',...
                    'FunctionTolerance',1e-6,...
                    'MaxIterations',2000,...
                    'SearchFcn','MADSPositiveBasis2N',...
                    'UseParallel',true,...
                    'PlotFcn',{'psplotbestf','psplotbestx'});
        
%         optimproblem.options = optimoptions('surrogateopt',...
%             'MaxTime',300,...
%             'Display','iter',...
%             'UseParallel',true,...
%             'PlotFcn',{'surrogateoptplot'});
        
        [xval(k,:),fval(k)] = patternsearch(optimproblem);
    end
    
    OptimizationResults.(strrep(KO{a},'/','_')).xval = xval;
    OptimizationResults.(strrep(KO{a},'/','_')).fval = fval;
end

save(['Data/PS_All_combined_new_v2_',num2str(a),'.mat']);

%% Initiation
close all;clc;clear;
addpath('AUX Functions','Main Functions','Data');
load Data.mat
load GenericNetwork.mat

%% Step 3a. Define fitting parameters for Pattern Search Algorithm (stochastic global optimization)
% Please refer to the function descriptions of patternsearch and optimoptions 
% for detailed explanations. 

%% Step 3a. Define fitting parameters for Pattern Search Algorithm (stochastic global optimization)
% Please refer to the function descriptions of patternsearch and optimoptions 
% for detailed explanations. 

% Specify the names of glycoprofiles to be fitted
ProfSel  = {'WT'}; % a cell of strings of profile names selected to be fitted.
% The name must be present in GenericNetwork.ProfNames.

num = 2; % Number of models fitted for each profile
method = 'StepProgression'; % the method used to compute the stationary state of Markov models
[~,transport_idx] = ismember({'cg2mg','mg2tg','tg2ab'},GenericNetwork.RxnTypes); % indices of transports in RxnTypes
optimproblem.x0 = 10*rand(1,length(GenericNetwork.Geneidx))-5;optimproblem.x0(transport_idx) = 0; % initial point a random vector with each element in the range of [-5,5]; set transport reactions to 0
        optimproblem.lb = ones(1,length(GenericNetwork.Geneidx)).*-5;optimproblem.lb(transport_idx) = 0; % lower bound for each parameter; set transport reactions to 0
        optimproblem.ub = ones(1,length(GenericNetwork.Geneidx)).*5;optimproblem.ub(transport_idx) = 0; % upper bound for each parameter; set transport reactions to 0
        optimproblem.solver = 'patternsearch'; % choose pattern search algorithm
        
        % No other constraints
        optimproblem.Aineq = [];
        optimproblem.bineq = [];
        optimproblem.Aeq = [];
        optimproblem.beq = [];
        optimproblem.nonlcon = [];
        optimproblem.intcon = [];
        optimproblem.rngstate = [];
        
%% Step 3b. Fit Markov models to glycoprofiles by Pattern Search Algorithm (stochastic glob
load GenericNetwork.mat

%% Step 3a. Define fitting parameters for Pattern Search Algorithm (stochastic global optimization)
% Please refer to the function descriptions of patternsearch and optimoptions 
% for detailed explanations. 

% Specify the names of glycoprofiles to be fitted
ProfSel  = {'WT'}; % a cell of strings of profile names selected to be fitted.
% The name must be present in GenericNetwork.ProfNames.

num = 2; % Number of models fitted for each profile
method = 'StepProgression'; % the method used to compute the stationary state of Markov models
[~,transport_idx] = ismember({'cg2mg','mg2tg','tg2ab'},GenericNetwork.RxnTypes); % indices of transports in RxnTypes
optimproblem.x0 = 10*rand(1,length(GenericNetwork.Geneidx))-5;optimproblem.x0(transport_idx) = 0; % initial point a random vector with each element in the range of [-5,5]; set transport reactions to 0
        optimproblem.lb = ones(1,length(GenericNetwork.Geneidx)).*-5;optimproblem.lb(transport_idx) = 0; % lower bound for each parameter; set transport reactions to 0
        optimproblem.ub = ones(1,length(GenericNetwork.Geneidx)).*5;optimproblem.ub(transport_idx) = 0; % upper bound for each parameter; set transport reactions to 0
        optimproblem.solver = 'patternsearch'; % choose pattern search algorithm
        
        % No other constraints
        optimproblem.Aineq = [];
        optimproblem.bineq = [];
        optimproblem.Aeq = [];load GenericNetwork.mat

%% Step 3a. Define fitting parameters for Pattern Search Algorithm (stochastic global optimization)
% Please refer to the function descriptions of patternsearch and optimoptions 
% for detailed explanations. 

% Specify the names of glycoprofiles to be fitted
ProfSel  = {'WT'}; % a cell of strings of profile names selected to be fitted.
% The name must be present in GenericNetwork.ProfNames.

num = 2; % Number of models fitted for each profile
method = 'StepProgression'; % the method used to compute the stationary state of Markov models
[~,transport_idx] = ismember({'cg2mg','mg2tg','tg2ab'},GenericNetwork.RxnTypes); % indices of transports in RxnTypes
optimproblem.x0 = 10*rand(1,length(GenericNetwork.Geneidx))-5;optimproblem.x0(transport_idx) = 0; % initial point a random vector with each element in the range of [-5,5]; set transport reactions to 0
        optimproblem.lb = ones(1,length(GenericNetwork.Geneidx)).*-5;optimproblem.lb(transport_idx) = 0; % lower bound for each parameter; set transport reactions to 0
        optimproblem.ub = ones(1,length(GenericNetwork.Geneidx)).*5;optimproblem.ub(transport_idx) = 0; % upper bound for each parameter; set transport reactions to 0
        optimproblem.solver = 'patternsearch'; % choose pattern search algorithm
        
        % No other constraints
        optimproblem.Aineq = [];
        optimproblem.bineq = [];
        optimproblem.Aeq = [];