%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data','Data/OptimizationResults');
load Data.mat;
load GenericNetwork.mat;
load ProcessedModels.mat

% Select the profiles to visualize
ProfSel = fieldnames(OptimizationResults);

%% Step 5a: Visualize other profiles in comparison with the WT profile
% Visualize the performance of fitted models and compare them with random
% models

% VisualizePerformance(OptimizationResults,GenericNetwork);

% Input:
% 1. OptimizationResults: the struct variable containing the info
% regarding the fitted models and computed from Step 3.

VisualizePerformance(OptimizationResults);


 %% Step 5b. Compare transition probabilities,  pseudo fluxes, and TP sensitivities between models of two different glyprofiles

 % select the name of the base glycoprofile against which all other glycoprofiles are
 % compared, usually a wildtype glycoprofile.
ComparisonProf1 = 'WT';

for a = 1:length(ProfSel)
    if ~strcmp(ComparisonProf1,ProfSel{a})
        % Compare and visualize significantly perturbed reactions in terms of TP fold
        % changes, model pseudo-flux fold changes, and sensitivities of the
        % corresponding TPs. Significant differentiations were assigned yellow circles.
        % Significantly perturbed reaction types are highlighted red.

        % [TPComp,PseudoFluxComp,SensitivityComp,RxnComp] = CompareTPandPseudoFlux(Prof, BaseProf,OptimizationResults,GenericNetwork);

        % Input:
        % 1. Prof: the string of the selected profile name for
        % comparison
        % 2. BaseProf: the string of the selected of the base profile used for comparison
        % 3. OptimizationResults: the struct variable containing the info
        % regarding the fitted models and computed from Step 3.
        % 4. GenericNetwork: the struct variable containing the info regarding
        % a generic N-glycosylation network and computed from Step 2.

        % Output:
        % 1. TPComp: a struct containing the TP log10(fold change) (double-type variable) and
        % significance (double-type variable, ranksum test between all TP fold changes and sample permutations).
        % 2. PseudoFluxComp: a struct containing the median model pseudo-flux log10(fold change) (double-type variable) and
        % significance (double-type variable, ranksum test between all pseudo-flux fold changes and sample permutations).
        % 3. SensitivityComp: a struct containing the maximum sensitivity metrics computed from Step 4 (abs(RMSE%/Perturbation%))
        % and significance for each considered reaction types and compared profiles.
        % 4. RxnComp: a cell of strings containing the row labels of the
        % reaction types considered for TPComp, PseudoFluxComp, and SensitivityComp

        CompName = [ProfSel{a} '_vs_' ComparisonProf1];
        [ComparativeResults.(CompName).TPComp,...
            ComparativeResults.(CompName).PseudoFluxComp,...
            ComparativeResults.(CompName).SensitivityComp,...
            ComparativeResults.(CompName).RxnComp] = ...
            CompareTPandPseudoFlux(ProfSel{a}, ComparisonProf1,OptimizationResults,GenericNetwork);
    end
end