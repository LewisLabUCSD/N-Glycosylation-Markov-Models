%% Initiation
close all;clc;clear;
addpath('AUX Functions','Main Functions','Data','Data/OptimizationResults');
load Data.mat;
load GenericNetwork.mat;
OptimizationResults = LoadOptimizationResults('OptimizationResults_Steric'); % load all optimization results and combine them together

% Select the profiles to visualize
ProfSel = {'WT'};

% Visualize fitted model results for each selected profiles, sequentially
for a = 1:length(ProfSel)
    %% Step 4a. Compute model characteristics
    
    %%%%%%%%%%%%%%%%%%%% Computed fitted models %%%%%%%%%%%%%%%%%%%%
    % ComputeFittedModels simulate the Markov models with fitted parameters
    % whereas GenerateRandomModels simulate the Markov models with random
    % parameters.
    % OptimizationResults = ComputeFittedModels(Prof,GenericNetwork,DataSet,OptimizationResults);
    % RandomResults = GenerateRandomModels(ProfSel{a},GenericNetwork,DataSet,OptimizationResults,simNum);
    
    % Input:
    % 1. Prof: the string of the selected profile name for
    % computing fitted models
    % 2. GenericNetwork: the struct variable containing the info regarding
    % a generic N-glycosylation network and computed from Step 2
    % 3. DataSet: the struct variable containing the info regarding the
    % experimental data and computed from Step 1
    % 4. OptimizationResults: the struct variable containing the info
    % regarding the fitted models and computed from Step 3.
    
    % Output:
    % 1. OptimizationResults: struct variable constructed from Step 3 and
    %    contains information of fitted models. The following fields were
    %    appended to each profile:
    %    a. error: the error computed by the objective function
    %    b. Predata_noRes: predicted signal corresponding to m/z values in mz_all
    %    c. Predata_raw: predicted signal corresponding to the absorbed glycan isoforms in Glys
    %    d. Glys_raw: the absorbed glycan isoforms in Glys
    %    e. PseudoConc: pseudo-concentrations of the intermediate glycans
    %    f. Glys_conc: the glycan isoforms corresponding to PseudoConc entries
    %    g. PseudoFlux: pseudo-fluxes of network reactions
    %    h. Rxn_flux: the edge list (reactant, product) representing the network
    %    i. AnnotatedMz: m/z values with annotations
    %    j. AnnotedGlycans: glycan annotations at the annotated m/z values
    
    OptimizationResults = ComputeFittedModels(ProfSel{a},GenericNetwork,DataSet,OptimizationResults);
    RandomResults = GenerateRandomModels(ProfSel{a},GenericNetwork,DataSet,OptimizationResults,20);
    OptimizationResults.(ProfSel{a}).RandomResults = RandomResults; % record data
    
    %% Step 4b. visualize fitted transition probabilities
    % [xval_media,xval]= PlotTPbyComp(ProfSel{a},OptimizationResults);
    
    % Input:
    % 1. ProfSel{a}: the string of the selected profile name for
    % visualizing transition probabilities
    % 2. OptimizationResults: struct variable constructed from Step 3 and
    % contains information of fitted models. Please refer to Step 3 for detailed info.
    % 3. GenericNetwork: the struct variable containing the info regarding
    % a generic N-glycosylation network and computed from Step 2
    
    % Output:
    % 1. xval_media: the median log10(transition probability) for each reaction
    % type among all fitted models of a profile
    % 2. xval (not shown): the log10(transition probability) of all fitted models for each reaction
    % type of a profile
    [OptimizationResults.(ProfSel{a}).xval_median]= PlotTPbyComp(ProfSel{a},OptimizationResults,GenericNetwork);
    
    %% Step 4c. visualize fitted glycoprofiles
    % [plotData,plotErr] = PlotPredVsExp(ProfSel{a},OptimizationResults);
    
    % Input:
    % 1. ProfSel{a}: the string of the selected profile name for
    % visualizing transition probabilities
    % 2. OptimizationResults: struct variable constructed from Step 3 and
    %    contains information of fitted models. Please refer to Step 3 for detailed info.
    % 3. numSel (optional): number of top signals to be plotted. m/z with the
    % highest experimental signal intensities will be plotted. If not
    % specified, top 20 signals (or all signals if total m/z values with non-zero
    % signals) will be plotted.
    
    % Output:
    % 1. ExpVsPredData: the signal intensities for each selected signals to be plotted (m/z).
    % Each row represents a fitted model and each column represents a
    % specific m/z.
    
    OptimizationResults.(ProfSel{a}).ExpVsPredData = PlotPredVsExp(ProfSel{a},OptimizationResults,20);
    
    %% Step 4d. visualize glycoform ratios for each m/z
    % GlycoformData = PlotGlycoforms(ProfSel{a},OptimizationResults,GenericNetwork,20, 0.001);
    
    % Input:
    % 1. ProfSel{a}: the string of the selected profile name for
    % visualizing transition probabilities
    % 2. OptimizationResults: struct variable constructed from Step 3 and
    %    contains information of fitted models. Please refer to Step 3 for detailed info.
    % 3. GenericNetwork: the struct variable containing the info regarding
    % a generic N-glycosylation network and computed from Step 2.
    % 4. numSel (optional): number of top signals to be plotted. m/z with the
    % highest experimental signal intensities will be plotted. If not
    % specified, top 20 signals (or all signals if total m/z values with non-zero
    % signals) will be plotted.
    % 5. threshold: signal intesity threshold below which glycoforms are
    % ignored
    
    % Output:
    % 1. GlycoformData: glycoform data with the following field:
    %    a. plotData: the dataset used to plot the heatmap showing the
    %    relative glycoform ratios at each m/z
    %    b. mz: m/z values considered for the analysis
    %    c. AllAsGlys: all glycoforms considered for the analysis
    
    OptimizationResults.(ProfSel{a}).GlycoformData = PlotGlycoforms(ProfSel{a},OptimizationResults,GenericNetwork,20, 0.001);
    
    %% Step 4e. visualize model pseudo-fluxes of fitted glycoprofiles
    % [plotData,plotErr] = PlotPredVsExp(ProfSel{a},OptimizationResults);
    
    % Input:
    % 1. ProfSel{a}: the string of the selected profile name for
    % visualizing transition probabilities
    % 2. OptimizationResults: struct variable constructed from Step 3 and
    %    contains information of fitted models. Please refer to Step 3 for detailed info.
    % 3. GenericNetwork: the struct variable containing the info regarding
    % a generic N-glycosylation network and computed from Step 2.
    % 4. numSel (optional): number of top produced glycoforms (in terms of
    % pseudo-fluxes feeding into their respective absorption state)
    % to be considered for plotting the graph of major network fluxes. If not
    % specified, top 20 glycoforms will be considered (or all signals if total m/z values with non-zero
    % signals) will be plotted.
    
    % Output:
    % 1. FluxesbyComp: total model pseudo-fluxes for reaction types, in the
    % order specified by the variable RxnTypes.
    % 2. Subnetwork: struct variable containing the following info
    % regarding the visualized model through the synthesis network; has the following
    % fields:
    %   a. G_trimmed: graph object of the trimmed network used to visualize
    %   the major model fluxes leading to the syntheses of selected
    %   glycans.
    %   b. NodeWT/EdgeWt: pseudo-concentration/fluxes through the trimmed
    %   model, in the orders specified by G_trimmed.Nodes.Name and
    %   G_trimmed.Edges.EndNodes, respectively.
    
    [OptimizationResults.(ProfSel{a}).FluxesbyComp,OptimizationResults.(ProfSel{a}).Subnetwork] = PlotModelFluxes(ProfSel{a},OptimizationResults,GenericNetwork,20);
    
    %% Sensitivity analysis (slow, can take >10 min per analysis, comment out the function if do not need the analyses)
    % [plotData,plotErr] = PlotPredVsExp(ProfSel{a},OptimizationResults);
    
    % Input:
    % 1. ProfSel{a}: the string of the selected profile name for
    % visualizing transition probabilities
    % 2. GenericNetwork: the struct variable containing the info regarding
    % a generic N-glycosylation network and computed from Step 2.
    % 3. DataSet: the struct variable containing the info regarding the
    % experimental data and computed from Step 1
    % 4. OptimizationResults: struct variable constructed from Step 3 and
    %    contains information of fitted models. Please refer to Step 3 for detailed info.
    % 5. UseNumofSamples (optional): a scalar. If the sample size is too
    % big, users may specify a subset of samples to conduct the
    % perturbation analysis and save time. If UseNumofSamples is bigger 
    % than the total sample size, the total sample size will be used.
    
    % Output:
    % 1. SensitivityAnalysis: a struct variable containing the data of
    % sensitivity analysis. The variable has the following field:
    %   a. perturbationPct: percentage perturbation (+ or -) of a reaction's 
    % transition probabilities
    %   b. RxnNames: reaction types purturbed in the sensitivity analysis
    %   c. errorMat: the sensivitity values are computed by perturbed model
    %   RMSE divided by percentage perturbation. Each row represents a
    %   reaction type and each column represents a perturbation percentage.
    
    OptimizationResults.(ProfSel{a}).SensitivityAnalysis = ConductSensitivityAnalysis(ProfSel{a},GenericNetwork,DataSet,OptimizationResults,5);
    
end

save('Data/ProcessedModels.mat','OptimizationResults');