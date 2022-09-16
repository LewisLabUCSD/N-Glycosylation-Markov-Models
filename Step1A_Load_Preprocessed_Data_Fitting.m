%% Initiation
close all;clc;clear;
addpath('Aux Functions','Main Functions','Data');

DatafileName = 'Data.xlsx';

%% Step 1a. Load preprocessed data
%%%%%%%%%%%%%%%%%%%% Required Data for fitting %%%%%%%%%%%%%%%%%%%%%%%%%
[~,ProfNames] = xlsread(DatafileName,'MS Raw','1:1');ProfNames = ProfNames(3:end)'; ProfNames = strrep(ProfNames,'/','_');% name of each glycoprofile
mz_all = xlsread(DatafileName,'MS Raw','A:A'); % all measured m/z values

[~,compositions] = xlsread(DatafileName,'MS Raw','B:B'); compositions = compositions(2:end);% the monossacharide compositions
profiles = xlsread(DatafileName,'MS Raw');profiles = profiles(:,3:end); % relative intensities measured at the corresponding m/z
profiles(isnan(profiles)) = 0; % fill empty elements
% values, where each column represent a glycoprofile.
profiles = profiles./sum(profiles); % normalize the  sum of signal intensities in each profile to 1 

%%%%%%%%%%%%%%%%%%%% optional data for fitting: glycan linkage annotation %%%%%%%%%%%%%%%%%%%%%%%%%
mz = xlsread(DatafileName,'Annotation','A:A'); % annotated m/z values
[~,LinkageResStruct] = xlsread(DatafileName,'Annotation','B:B');LinkageResStruct = LinkageResStruct(2:end); % all annotated glycan structures from all profiles in the linear code format
LinkageResStructSel = xlsread(DatafileName,'Annotation');LinkageResStructSel = logical(LinkageResStructSel(:,3:end)); % selection flag of annotated glycan structures for each profile 
LinkageInfoAvail = any(LinkageResStructSel~=0); % whether glycan linkage annotation is present or not for each profile
DataSet = ws2struct;

%%%%%%%%%%%%%%%%%%%% visualize experimental data for sanity check %%%%%%%%%%%%%%%%%%%%%%%%%

% visualizeExpData(DataSet, ProfName, signalSel)

% Input:
% 1. DataSet: DataSet variable created from the previous steps in the
% script
% 2. ProfName: a string or a cell vector of strings of glycoprofile name(s)
% to be plotted. ProfName should be present in the variable ProfNames. If
% users would like to plot all glycoprofiles, enter [].
% 3. signalSel: specify the number of top glycan structure annotations to be plotted.
% enter [] to plot all annotations.

% Output:
% 1. Figure type 1: relative glycomic signals (glycoprofiles)
% 2. Figure type 2: glycan structures at annotated m/z
visualizeExpData(DataSet,[],[]);

save('Data/Data.mat','DataSet');
