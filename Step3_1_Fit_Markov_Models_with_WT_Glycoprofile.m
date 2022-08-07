%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data');
load Data.mat
load GenericNetwork.mat

%% Step 3a. Define fitting parameters for Pattern Search Algorithm (stochastic global optimization)
% Please refer to the function descriptions of patternsearch and optimoptions
% for detailed explanations. Specify the names of glycoprofiles to be
% fitted. The names provided must be present in GenericNetwork.ProfNames.

ProfSel  = {'WT'}; % a cell of strings of profile names selected to be fitted.

%% Step 3b. Fit Markov models to glycoprofiles by stochastic global optimization

% initiate a struct to store fitted parameters
OptimizationResults = struct;

% Each selected glycoprofiles is fitted sequentially
for a = 1:length(ProfSel)

    num = 2; % Number of models fitted for each profile

    for k = 1:num

        %%%%%%%%%%%%%%%%%%%%% Set up optimization problem %%%%%%%%%%%%%%%%%%%%%
        % For advanced users, please open the function script SetUpFittingProblem
        % to review detailed set-up instruction and modify fitting parameters.
        % Default setup is used here. The default setup was used in our
        % previous work.

        % OptimizationProblem = SetUpFittingProblem(Num,GenericNetwork,DataSet);

        % Input:
        % 1. Prof: the string of the selected profile name for fitting (must also be an element in DataSet.ProfNames
        % 2. DataSet: struct variable contructed from Step 1.
        % 3. GenericNetwork: struct variable constructed from Step 2.
        % 4. StericFlag: true or false, whether to consider the impact of
        % steric interaction. Default is false.
        % 5. UseWTStericFlag: whether to apply the fitted WT steric factors to
        % other profiles. StericFlag must be true if this parameter is set
        % true. If false, the algorithm will fit a  set of steric factors
        % (SFs). Default is false.
        % for the profile
        % 6. WTSteric: a double-type vector specifying the the wild type steric factors. If WTStericFlag
        % is false, set WTSteric as []. In the current default setting, the SFs will be
        % fitted if fitting a wildtype profile, whereas the wildtype SFs will
        % be fed to the fitting algorithm if fitting a glycoengineered profile.
        % 7. Method: algorithm chosen for stochastic global optimizations.
        % Choose beteween 'PatternSearch' or 'SAPSHybrid'.

        % Output:
        % 1. Optimization Problem: a struct variable with the following field:
        %    a. optimproblem: actual optimproblem structure containing the fed
        %    to the function "patternsearch". Refer to MathWork for more info.
        %    b. mzRes: double-type indices of m/z values at which glycan structures were annotated
        %    c. linkagePosRes: cell of strings representing the linear codes of actual glycan annotations
        %    d. NumOfModels: numbers of models fitted for the selected profiles
        %    e. ExpData: experimentally measured glycoprofiles (normalized)
        %    f. xval: pre-initiated variable storing the fitted transition probabilities for each "reactions" in RxnTypes
        %    g. fval: pre-initiated variable storing the fitting errors


        StericFlag = false;
        UseWTStericFlag = false;
        WTSteric = [];
        Method = 'PatternSearch';

        OptimizationProblem = SetUpFittingProblem(num,ProfSel{a},GenericNetwork,DataSet,StericFlag,UseWTStericFlag,WTSteric,Method);


        %%%%%%%%%%%%%%%%%%%%% Fitting (global optimization) %%%%%%%%%%%%%%%%%%%%%
        % OptimizationResults = runGlobalOptimization(Prof, OptimizationResults,OptimizationProblem);

        % Input:
        % 1. Prof: the string of the selected profile name for fitting (must also be an element in DataSet.ProfNames
        % 2. DataSet: struct variable contructed from Step 1.
        % 3. OptimizationProblem: struct variable constructed from
        % SetUpFittingProblem.

        % Output:
        % 1. OptimizationResults: struct variable storing the optimized
        % parameters and errors

        OptimizationResults = runGlobalOptimization(ProfSel{a}, OptimizationResults,OptimizationProblem);
    end
end

%% Step 3c. Store the fitting result
save('Data/OptimizationResults/OptimizationResults_WT_42.mat','OptimizationResults');