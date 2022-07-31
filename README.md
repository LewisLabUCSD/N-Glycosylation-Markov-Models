
# **Toward rational glycoengineering: Markov models for *N*-Glycosylation**

## **Introduction**
Glycosylated biopharmaceuticals are important in the global pharmaceutical market. Despite the importance of their glycan structures, our limited knowledge of the glycosylation machinery still hinders the controllability of this critical quality attribute. To facilitate discovery of glycosyltransferase specificity and predict glycoengineering efforts, here we extend the approach to model N-linked protein glycosylation as a Markov process. Our model leverages putative glycosyltransferase (GT) specificity to define the biosynthetic pathways for all measured glycans, and the Markov chain modeling is used to learn glycosyltransferase isoform activities and predict glycosylation following glycosyltransferase knock-in/knockout.

## **Requirement**
The pipeline is based on our previously published [work](https://www.sciencedirect.com/science/article/pii/S2590262820300010). **MATLAB** **ver. 2019b** or above is required to ensure proper execution of all functions. Operators and Elementary Operations, Global Optimization, Parallel Computing, Statistics and Machine Learning, Econometrics, Curve Fitting, and System Identification, Graph and Network Algorithms Toolboxes. All third-party functions were cited as comments at the beginnings of corresponding code scripts. Minimal **16 GB RAM** and two **4-core** CPU(s) are highly recommended to ensure efficient running and fitting processes. 

## **Pipeline Overview**
The basic pipeline of model fitting and visualization include the following steps:
- **Step 0: Data preparation**
  - Raw glycomics data and annotations (if available) need to be gathered and processed prior to fitting. Glycomics data with glycan structures/glycan compositions but no m/z values can also be used. If you do not have the model-compatible glycan compositions but only the glycan structures, ***Sup1_Get_compositions_from_linearcodes.m*** in the ***Supplemental Steps*** folder can be used to obtain corresponding glycan compositions from user-provided linear codes. Annotations of high-intensity or branching glycoforms are highly recommended.
  - ***Data.xlsx*** in the Data folder is the only input required to run the basic pipeline (**Figure 1**). Please refer to ***Step0_Data_Preprocessing.pptx*** for details on how to organize your data.
  - Example glycoprofiles (used in our previous [work](https://www.sciencedirect.com/science/article/pii/S2590262820300010)) were included in the example  ***Data.xlsx*** 

- **Step 1: Load Preprocessed Data**
  - This step reads the data from ***Data.xlsx*** and stores the dataset as a structure variable *DataSet*, which is used for fitting & analyses. 
  - Directly run the script ***Step1_Load_Preprocessed_Data.m*** to obtain ***Data.mat***, a MATLAB file containing *DataSet* and stored in ***Data*** folder. The existing *Data.mat* was generated for the example glycoprofiles used for fitting demonstration.
  - For practicality, the complexity level of the network should not exceed 23.
 
 - **Step 2: Generate a generic *N*-glycosylation Network**
   - This step generates a *N*-glcosylation synthetic network by recursively modifying the high-mannose structure (Man9) using reaction rules as specified in our previous [work](https://www.sciencedirect.com/science/article/pii/S2590262820300010). The size of the network can be controlled based on the scope of glycan structures considered or an arbitrary complexity level. The network is organized as an adjacency matrix for constructing Markov models.
   - Directly run the script ***Step2_Create_Genric_Glycan_Network.m*** to obtain ***GenericNetwork.mat***, a MATLAB file containing *GenericNetwork* and stored in ***Data*** folder. *GenericNetwork* contains all information regarding the generic model used for fitting.

 - **Step 3_1/2: Fit a generic Markov model to WT/other glycoprofiles**
   - This step generates fitted  Markov models for each specified glycoprofile. Using the default setting, the time required to obtain models for all the glycoprofile is approximately 3 hr/optimization x number of models fitted for each profile x number of profiles/total MATLAB sessions. Please make sure that CPU/memory utilization is not consistently at 100% when the optimization is running. Please refer to MathWork for more details on parallel processing. 
   - Fitting the WT glycoprofile first is required when the model steric parameters (Change Log Jul-26-2022) are considered for model fitting (running ***Step3_1_Fit_Markov_Models_with_WT_Glycoprofile.m*** first and ***Step3_2_Fit_Markov_Models_with_other_Glycoprofiles.m*** second). If steric parameters are not considered, either scripts can be used. To obtain more reliable steric factors, users are suggested to fit as many wildtype models as pratically possible.
   - Before running either script, specify the variables *ProfSel* (model names selected to be fitted),  *StericFlag* (true or false, whether to consider the steric impact), and *UseWTStericFlag* (true or false, whether to use the steric impact learnt from the wild type models). Refer to comments in the script for more details.
   - For debug, the evaluation of the Pattern Search Algorithm is plotted in real time and printed to the command window. (**Figure 1**). 
   - Upon completion, running either script will start the fitting process and store the fitted model parameters in *OptimizationResults*. For consistent organization when running multiple MATLAB sessions, the *OptimizationResults* from each session will be stored in ***yourPreferredName_#d.mat*** files (where #d is a number) in the folder ***Data/OptimizationResults/***. Only one *yourPreferredName* should be used for all glycoprofiles in the same study to ensure our customized file loaders can load all relevant files. ***OptimizationResults_Steric_WT_#d.mat*** and ***OptimizationResults_Steric_others_#d.mat*** were generated with the example dataset.

![enter image description here](https://github.com/LewisLabUCSD/N-Glycosylation-Markov-Models/blob/main/Figures/Figure%201.PNG)

**Figure 1**: (Top) Values of objective function at each optimization iteration with Pattern Search Algorithm. (Bottom) The currently best set of model parameters (transition probabilities associated with reaction types) that minimizes the objective function. 

 - **Step 4: Visualize predicted glycoprofiles from fitted Markov models**
   - ***Step4_Visualize_Fitted_Glycoprofiles.m*** quantifies a variety of model features (e.g. pseudo-fluxes, pseudo-concentration) and produces visuals by simulating the generic models with the fitted parameters from Step 3. *ProfSel* should be specified first to selected the glycoprofiles for visualization. A few example visuals are shown in **Figure 2** & **Figure 3**.   
 

![enter image description here](https://github.com/LewisLabUCSD/N-Glycosylation-Markov-Models/blob/main/Figures/Figure%202.PNG)

**Figure 2**: Visualization of model features from models fitted with the wildtype *N*-glycoprofile. (A) The model fluxes through each of the reaction types. (B) The transition probabilities (model parameters) associated with each of the reaction types, including the steric factors. (C) Comparison between the experimental glycoprofile and the predicted glycoprofiles produced from the fitted Markov models in terms of realtive intensities (MS signals). (D) The relative ratios (shades of color) of major glycoforms at each m/z values predicted by the fitted models. Red glycans means it is the most abundant glycoform at a specific m/z value.   

![enter image description here](https://github.com/LewisLabUCSD/N-Glycosylation-Markov-Models/blob/main/Figures/Figure%203.PNG)

**Figure 3** Visualization of the network topology of the from models fitted with the wildtype *N*-glycoprofile. Each node represents a glycan and each edge represents a glycosyltransferase/glycosidase reaction (reaction types). (A) The active reactions and intermediate glycans (blue edges and blue nodes) leading to the production of major glycoforms (red nodes), in the context of a generic network (grey edges and nodes). (B) The major fluxes (red edges) leading to the productiong of major glycoforms (red nodes), in the context of the active network (as show in (A)). 

### Please refer to the comments in the MATLAB scripts for more details on:
- Individual steps (functions) and the variables (including their formats) required to run the pipelines
- Customizing fitting parameters or visualization
- Accessing fitted model parameters and features for additional analyses

## **Change Log**

- Jul-26-2022
  - Updated consturction, fitting, and visualization functions to make them compatible when users choose to consider the impact of steric interactions of *N*-glycan antenna on model parameters (transition probabilities):
    - Rxns considered included a3SiaT, b4GalT, iGnT, GnTIV, and GnTV.
    - Quantified impact of steric interactions to each considered rxn as a parameter.    

- Jul-20-2022 
  - Removed constraint of fucosylation (a6FucT) to include all theoretically possible glycans (*CreateGlycanRxnList.mat*). 
  - Restricted the maximal length of LacNAc chain in the function *CreateGlycanRxnList*.

## **Reference**
1. [A Markov model of glycosylation elucidates isozyme specificity and glycosyltransferase interactions for glycoengineering](https://www.sciencedirect.com/science/article/pii/S2590262820300010\).)

## **Additional Info**
For additional questions, please reach out to [Systems Biology and Cell Engineering – Lewis Lab](https://lewislab.ucsd.edu/) at the University of California, San Diego.  
