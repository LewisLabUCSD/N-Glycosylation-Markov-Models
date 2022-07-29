%% Initiate environment
addpath('../AUX Functions','../Main Functions');

%% Copy and paste your linearcodes into the variable linearcodes

% An example of linearcodes variable is provided here. To specify your own
% linear codes, the easiest way is to create an empty cell named
% linearcodes and open this variable from the Workspace. Copy and paste
% your linearcode data from your own dataset (best from Excel) directly into
% the opened variable.

linearcodes = {'(Ab4GNb2Ma3(Ab4GNb2Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(Ab4GNb2Ma3(NNa3Ab4GNb2Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(NNa3Ab4GNb2Ma3(NNa3Ab4GNb2Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa3Ab4GNb2Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(Ab4GNb2(Ab4GNb4)Ma3(NNa3Ab4GNb2(NNa3Ab4GNb6)Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa3Ab4GNb2(NNa3Ab4GNb6)Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(NNa3Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa3Ab4GNb2(NNa3Ab4GNb6)Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(NNa3Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa3Ab4GNb2(NNa3Ab4GNb3Ab4GNb6)Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(NNa3Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa3Ab4GNb2(NNa3Ab4GNb3Ab4GNb3Ab4GNb6)Ma6)Mb4GNb4(Fa6)GN);Asn';
    '(NNa3Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa3Ab4GNb2Ma6)Mb4GNb4(Fa6)GN);Asn'};

%% Get Glycan compositions

% compositions contains the pipeline-compatible notations of glycan
% compositions. You may open the variable in the Workspace and copy and
% paste the composition data back to your dataset.
compositions = GetGlycanCompositions(linearcodes);
