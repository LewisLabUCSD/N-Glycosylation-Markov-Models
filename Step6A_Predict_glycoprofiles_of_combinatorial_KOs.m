%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data','Data/OptimizationResults');
load Data.mat;
load GenericNetwork.mat;
load ProcessedModels.mat;
load ComparativeResults.mat

%% Step 6a. Predict Combinatorial knockout from single KO fold changes
% Fitted widetype and single-knockout models are required to predict the
% combinatorial knockout glycoprofiles

KnockoutSel = {{'St3gal4','St3gal6'};...
    {'Mgat2','St3gal4','St3gal6'}};
BaseProfSel = 'WT';

for a = 1:length(KnockoutSel)

    ProfName = strjoin(KnockoutSel{a},'_');

    % Predict the model parameters for combinatorial knockout
    [PredictedResults.(ProfName).xval,PredictedResults.(ProfName).LacNAcLenPenalty] = ...
        PredictCombKOProfiles(KnockoutSel{a},BaseProfSel,ComparativeResults,OptimizationResults,GenericNetwork);

    PredictedResults = ComputePredictedModels(ProfName,GenericNetwork,DataSet,PredictedResults);

    [PredictedResults.(ProfName).xval_median]= PlotTPbyComp(ProfName,PredictedResults,GenericNetwork);

    threshold = 1e-3;
    OptimizationResults.(ProfName).GlycoformData = PlotGlycoforms(ProfName,PredictedResults,GenericNetwork,15, threshold);

    numSel = 20;
    OptimizationResults.(ProfName).ExpVsPredData = PlotPredGlycoprofile(ProfName,PredictedResults,numSel);

    % Comparing predicted results and the experimental results as a sanity
    % check, if the experimental results exist.
    % visualizeExpData(DataSet,ProfName,numSel);

end