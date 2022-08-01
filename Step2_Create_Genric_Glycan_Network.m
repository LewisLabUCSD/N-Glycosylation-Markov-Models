%% Initiation
close all;clc;clear;
addpath('AUX Functions','Main Functions','Data');
load Data.mat

%% Step 2a. Construct a generic N-glycosylation network;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[rxnList,Glys] = CreateGlycanRxnList(N,targetGlycans)

% Input:
% 1. N: complexity level, the number of total steps by
% glycosyltransferease/glycosidases to process a Man9 glycan. The network
% will include all intermediate glycans and reaction steps to be generated within N steps.
% 2. RxnSel: reaction types considered for constructing the biosynthetic
% N-glycosylation network. Specify as the string "default" to include all
% reactions.
% 3. targetGlycans(optional): a cell list of strings representing the
% linear codes of glycans which the user intend to include in the network.
% If all the target glycans are produced in the network at a complexity level < N,
% network extention will terminate at the completion of current complexity level.
% Network extention will terminate once N complexity level is reached
% regardless whether all target Glycans have been made by the network.

% Output:
% 1. AllrxnList: the edge list of all glycosylation reactions, intercompartal transport,
% and Glycan "secretion" included in the network. Note that "self absoprtions" of glycans,
% elements of Markov model, are not included in the list. Columns
% 1-6 are reactants (string), products (string), involved reaction types (string), involved compartment(s) (string),
% complexity levels (double) at which reactions happen, and unique reaction labels (string), respectively.
% 2. Glys: a cell vector of strings representing all intermediate glycans produced in the network, also
% distinguished by Golgi compartments and absorption state ("secreted").
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 23;
RxnSel = {'N/A' 'ManI' 'ManII' 'GnTI' 'GnTII' 'GnTIV' 'GnTV' 'a6FucT' 'b4GalT' 'a3SiaT' 'iGnT'};
[AllrxnList,Glys] = CreateGlycanRxnList(N,RxnSel);

%% Step 2b. Construct Transition Probability Matrix (TPM) from the generic N-glycosylation network

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[TM, AllrxnList_TMrow, AllrxnList_TMcol, AllrxnList_Br] = ConstructTPM(AllrxnList)

% Input:
% 1.AllrxnList: the edge list of all glycosylation reactions, intercompartal transport,
% and Glycan "secretion" included in the network. Note that "self absoprtions" of glycans,
% edge elements of Markov model, are not included in the list. Columns
% 1-6 are reactants (string), products (string), involved reaction types (string), involved compartment(s) (string),
% complexity levels (double) at which reactions happen, and unique reaction labels (string), respectively.
% 2. Glys: a cell vector of strings representing all intermediate glycans produced in the network
% specified by AllrxnList.

% Output:
% 1. TM: sparse double-type matrix, transition probability matrix (an adjacency matrix representing all 
% the edges from AllrxnList). Each row represents a reactant glycan and each
% column represents a product glycan. Each row is normalized to a sum of 1.
% The order of the row and column labels representing intermediate glycans are identical to the order of
% glycans in Glys. 
% 2. AllrxnList_TMidx: double-type vector; in the order of AllrxnList, the
% indices of elements in TM representing the reactions.  
% 3. AllrxnList_Br: double-type vector; in the order of AllrxnList, the
% branch of an N-glycan at which a monosshacharide is added. Only applicable to a3SiaT,
% b4GalT, and iGnT reaction types (otherwise assigned 0).  
% 4. AllrxnList_RxnTypes: cell of strings; in the order of AllrxnList, the
% antenna-specific (if applicable) reaction types associated with the reactions  
% 5. RxnTypes: cell vector of strings; all antenna-specific reaction types described by the model (Liang & Chiang et al., 2020)
% 6. AllrxnList_steric: whether a neighboring antenna on the reactant
% glycan is present (true or false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TM, AllrxnList_TMidx, AllrxnList_Br, AllrxnList_RxnTypes, RxnTypes,AllrxnList_steric,stericRxns] = ConstructTPM(AllrxnList,Glys);

%% Step 2c. Predefine generic variables for Markov model fitting with Pattern Search Algorithm

% Geneidx: cell vectors of double-type vectors; in the order of RxnTypes,
% each cell contains all the indices in TPM corresponding to the reactions of RxnTypes  
Geneidx = cellfun(@(x) AllrxnList_TMidx(strcmp(x,AllrxnList_RxnTypes)), RxnTypes,'UniformOutput',0);

% Predefine Pattern Search function constants
AbsGlys = Glys(contains(Glys,'[ab]')); % cell vector of strings of all absorbed glycans
NA = length(AbsGlys); % number of model nodes that are absorbed glycans
NT = length(Glys) - NA; % number of model nodes that are not absorbed glycans 
feedflux = 1; % total model flux feeding into the root glycan (Man9). 
pi0 = zeros(1,length(Glys)); % initial distribution state for all the glycans in the Markov model
pi0(strcmp(Glys,'(Ma2Ma2Ma3(Ma2Ma3(Ma2Ma6)Ma6)Mb4GNb4GN);Asn[cg]')) = feedflux; % assign Man9 (root glycan) with the total model flux
pi0_T = pi0(1:NT); % initial distribution state for all the glycans in the Markov model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AbsGlyIdx = GenerateAbsGlyIdx(NT, AbsGlys, compositions,mz_all)

% Input:
% 1. NT: the number of model nodes that are not absorbed glycans
% 2. AbsGlys: cell vector of strings of all absorbed glycans
% 3. compositions: cell of strings loaded from the data Excel file 
% representing the monossacharide compositions 
% 4. mz_all: double-type vector representing all the measured m/z values. Note
% that mz_all and composition are in the same corresponding order.

% Output:
% 1. AbsGlyIdx: a cell vector of double-type vectors; each cell represents a
% m/z value or a monossacharide composition. Multiple glycan isoforms are possible 
% for each m/z value or a monossacharide composition, and their position
% indices in the variable Glys are listed in the corresponding cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AbsGlyIdx = GenerateAbsGlyIdx(NT, AbsGlys, DataSet.compositions, DataSet.mz_all);

%% Save the generic model
GenericNetwork = ws2struct; 
save('Data\GenericNetwork_newNetwork.mat','GenericNetwork');
