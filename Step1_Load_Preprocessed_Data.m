%% Initiation
close all;clc;clear;
addpath('AUX Functions','Main Functions','Data');

%% Step 1a. Load preprocessed data
%%%%%%%%%%%%%%%%%%%% Required Data for fitting %%%%%%%%%%%%%%%%%%%%%%%%%
[~,ProfNames] = xlsread('Data.xlsx','MS Raw','1:1');ProfNames = ProfNames(3:end)'; % name of each glycoprofile
mz_all = xlsread('Data.xlsx','MS Raw','A:A'); % all measured m/z values
[~,compositions] = xlsread('Data.xlsx','MS Raw','B:B'); compositions = compositions(2:end);% the monossacharide compositions
profiles = xlsread('Data.xlsx','MS Raw');profiles = profiles(:,3:end); % relative intensities measured at the corresponding m/z
% values, where each column represent a glycoprofile.
profiles = profiles./sum(profiles); % normalize the  sum of signal intensities in each profile to 1 

%%%%%%%%%%%%%%%%%%%% optional data for fitting: glycan linkage annotation %%%%%%%%%%%%%%%%%%%%%%%%%
mz = xlsread('Data.xlsx','Annotation','A:A'); % annotated m/z values
[~,LinkageResStruct] = xlsread('Data.xlsx','Annotation','B:B');LinkageResStruct = LinkageResStruct(2:end); % all annotated glycan structures from all profiles in the linear code format
LinkageResStructSel = xlsread('Data.xlsx','Annotation');LinkageResStructSel = logical(LinkageResStructSel(:,3:end)); % selection flag of annotated glycan structures for each profile 
LinkageInfoAvail = any(LinkageResStructSel~=0); % whether glycan linkage annotation is present or not for each profile
DataSet = ws2struct;

save('Data/Data.mat','DataSet');
