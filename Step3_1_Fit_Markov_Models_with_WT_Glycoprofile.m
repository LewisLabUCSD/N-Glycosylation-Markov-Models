%% Initiation
close all;clc;clear;
addpath('AUX Functions','Main Functions','Data');
load Data.mat
load GenericNetwork_newNetwork.mat

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

    num = 4; % Number of models fitted for each profile

    for k = 1:num

        %%%%%%%%%%%%%%%%%%%%% Set up optimization problem %%%%%%%%%%%%%%%%%%%%%
        % For advanced users, please open the function script SetUpFittingProblem
        % to review detailed set-up instruction and modify fitting parameters.
        % Default setup is used here.

        % OptimizationProblem = SetUpFittingProblem(Num,GenericNetwork,DataSet);

        % Input:
        % 1. Prof: the string of the selected profile name for fitting (must also be an element in DataSet.ProfNames
        % 2. DataSet: struct variable contructed from Step 1.
        % 3. GenericNetwork: struct variable constructed from Step 2.
        % 4. StericFlag: true or false, whether to consider the impact of steric interaction
        % 5. UseWTStericFlag: whether to apply the fitted WT steric factors to
        % other profiles. StericFlag must be true if this parameter is set
        % true. If false, the algorithm will fit a new set of steric factors
        % (SFs)
        % for the profile
        % 6. WTSteric: a double-type vector specifying the the wild type steric factors. If WTStericFlag
        % is false, set WTSteric as []. In the current default setting, the SFs will be
        % fitted if fitting a wildtype profile, whereas the wildtype SFs will
        % be fed to the fitting algorithm if fitting a glycoengineered profile.
        %

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
        OptimizationProblem = SetUpFittingProblem(num,ProfSel{a},GenericNetwork,DataSet,StericFlag,UseWTStericFlag,[]);

        %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%
        repeatflag = true; % repeat the fitting if algorithm failed to converge in rare cases
        while repeatflag
            [xval,fval,~,output] = patternsearch(OptimizationProblem.optimproblem);
            if output.iterations > 1000
                repeatflag = false;
            end
        end

        %%%%%%%%%%%%%%%%%%%%% Record Fitting Results %%%%%%%%%%%%%%%%%%%%%

        % OptimizationResults = RecordOptimizationResults(ProfSel{a}, OptimizationResults, xval, fval,OptimizationProblem);

        % Input
        % 1. Prof: the string of the selected profile name for fitting
        % 2. OptimizationResults: struct variable used to store all
        % optimization results.
        % 3. xval: optimized model parameters of the current run
        % 4. fval: minimized objective function error of the current run
        % 5. OptimizationProblem: struct storing the fitting problem set up
        % created by the function SetUpFittingProblem.

        % Output
        % c. OptimizationResults: updated OptimizationResults variable. Has
        % the following fields:
        %   a.xval: double-type fitted transition probablities. Each row represent a fitted transition
        %     probabilities vectors and each column represents a reaction type (in the order of variable RxnTypes).
        %     Therefore, the number of rows is equal to the number of fitted
        %     models for a specific profile, whereas the number of columns is
        %     equal to the the length of variable RxnTypes.
        %   b.fval: the errors computed by the objective function for each
        %   fitted model.
        %   c. OptimizationProblem: the struct variable storing the fitting
        %   problem setups

        OptimizationResults = RecordOptimizationResults(ProfSel{a}, OptimizationResults, xval, fval,OptimizationProblem);
    end
end

%% Step 3c. Store the fitting result
save('Data/OptimizationResults/OptimizationResults_Steric_WT_NewNewNewNetWork_5.mat','OptimizationResults');