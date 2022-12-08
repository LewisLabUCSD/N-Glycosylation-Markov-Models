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

KnockoutSel = {{'Mgat2','St3gal4','St3gal6'}};
BaseProfSel = 'WT';

for a = 1:length(KnockoutSel)

    ProfName = strjoin(KnockoutSel{a},'_');

    % Predict the model parameters for combinatorial knockout
    [PredictedResults.(ProfName).xval,PredictedResults.(ProfName).LacNAcLenPenalty] = ...
        PredictCombKOProfiles(KnockoutSel{a},BaseProfSel,ComparativeResults,OptimizationResults,GenericNetwork);

    PredictedResults = ComputePredictedModels(ProfName,GenericNetwork,DataSet,PredictedResults);

    [PredictedResults.(ProfName).xval_median]= PlotTPbyComp(ProfName,PredictedResults,GenericNetwork);

    threshold = 1e-3;
    OptimizationResults.(ProfName).GlycoformData = PlotGlycoforms(ProfName,PredictedResults,GenericNetwork,25, threshold);

    numSel = 25;
    OptimizationResults.(ProfName).ExpVsPredData = PlotPredGlycoprofile(ProfName,PredictedResults,numSel);

end