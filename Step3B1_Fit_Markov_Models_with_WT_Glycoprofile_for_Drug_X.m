%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data');
load DrugXData.mat
load GenericNetwork_DrugX.mat

ProfSel  = {'WT'}; % a cell of strings of profile names selected to be fitted.

%% Step 3b. Fit Markov models to glycoprofiles by stochastic global optimization

% initiate a struct to store fitted parameters
OptimizationResults = struct;

% Each selected glycoprofiles is fitted sequentially
for a = 1:length(ProfSel)

    num = 5; % Number of models fitted for each profile

    for k = 1:num


        StericFlag = false;
        UseWTStericFlag = false;
        WTSteric = [];
        Method = 'PatternSearch';

        OptimizationProblem = SetUpFittingProblem(num,ProfSel{a},GenericNetwork,DataSet,StericFlag,UseWTStericFlag,WTSteric,Method);
        OptimizationResults = runGlobalOptimization(ProfSel{a}, OptimizationResults,OptimizationProblem);
    end
end

%% Step 3c. Store the fitting result
save('Data/OptimizationResults/OptimizationResults_DrugXWT_1.mat','OptimizationResults');