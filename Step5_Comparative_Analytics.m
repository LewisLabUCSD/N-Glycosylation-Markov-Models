%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data','Data/OptimizationResults');
load Data.mat;
load GenericNetwork.mat;
% load all optimization results and combine them together. the string input
% should be the file names excluding the last underscore and the number
% (_#d).
OptimizationResults = LoadOptimizationResults({'OptimizationResults_WT_','OptimizationResults_others_'}); 

% Select the profiles to visualize
ProfSel = fieldnames(OptimizationResults);

%% Visualize other profiles in comparison with the WT profile 

% Visualize fitted model results for each selected profiles, sequentially
for a = 1:length(ProfSel)
    
end