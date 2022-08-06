function [AbsGlyIdx,LeakageGlyIdx] = GenerateAbsGlyIdx(NT, AbsGlys, composition,mz_all)

% add up same monosaccharide types
generaltypes = {'Hex','HexNAc','NeuAc','dHex'}; % col

% calculate the total number of each general type of monosaccharide

composition_monoCounts_exp = zeros(length(composition),length(generaltypes));

for k = 1:length(composition)
    
    str = composition{k};
    
    % for dHex
    idx1 = strfind(str,'dHex');
    if ~isempty(idx1)
        flag = isstrprop(str(idx1+4),'digit');
        if ~flag
            composition_monoCounts_exp(k,4) = 1;
        else
            composition_monoCounts_exp(k,4) = str2double(str(idx1+4));
        end
    end
    
    % for NeuAc
    idx2 = strfind(str,'NeuAc');
    if ~isempty(idx2)
        flag = isstrprop(str(idx2+5),'digit');
        if ~flag
            composition_monoCounts_exp(k,3) = 1;
        else
            composition_monoCounts_exp(k,3) = str2double(str(idx2+5));
        end
    end
    
    % for HexNAc
    idx3 = strfind(str,'HexNAc');
    if ~isempty(idx3)
        flag = isstrprop(str(idx3+6),'digit');
        if ~flag
            composition_monoCounts_exp(k,2) = 1;
        else
            composition_monoCounts_exp(k,2) = str2double(str(idx3+6));
        end
    end
    
    % for Hex
    idx4 = strfind(str,'Hex');
    idx4 = setdiff(idx4,[idx1+1,idx3]);
    if ~isempty(idx4)
        flag = isstrprop(str(idx4+3),'digit');
        if ~flag
            composition_monoCounts_exp(k,1) = 1;
        else
            composition_monoCounts_exp(k,1) = str2double(str(idx4+3));
        end
    end
end

% Get general monosaccharide notations for Glys
composition_monoCounts_pre = zeros(length(AbsGlys),4);
for k1 = 1:length(AbsGlys)
    str = AbsGlys{k1};
    % Hex (col 1)
    composition_monoCounts_pre(k1,1) = length(strfind(str,'A'))-1 + length(strfind(str,'M'));% discount for Asn
    % HexNAc (col 2)
    composition_monoCounts_pre(k1,2) = length(strfind(str,'GN'));
    % NeuAc (col 2)
    composition_monoCounts_pre(k1,3) = length(strfind(str,'NN'));
    % dHex (col 2)
    composition_monoCounts_pre(k1,4) = length(strfind(str,'F'));
end
[~,Locb] = ismember(composition_monoCounts_pre,composition_monoCounts_exp,'rows');
AbsGlyIdx = arrayfun(@(x) find(Locb == x)+NT,(1:length(mz_all))','UniformOutput',0);

% Identify leakage compositions
Locd = ~Locb;
leakageComposition = unique(composition_monoCounts_pre(Locd,:),'rows');
LeakageGlyIdx = cell(size(leakageComposition,1),1);
for a = 1:size(leakageComposition,1)
    [~,Loc] = ismember(composition_monoCounts_pre,leakageComposition(a,:),'rows');
    LeakageGlyIdx{a} = find(Loc)+NT;
end

end