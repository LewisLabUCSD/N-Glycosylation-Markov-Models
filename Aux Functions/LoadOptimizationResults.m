function combResults = LoadOptimizationResults(str)

% Identify optimization result files from the Data folder
listing = dir('**/Data/OptimizationResults');
names = {listing.name};
names_all = {''};
for a = 1:length(str)
    names_temp = names(cellfun(@(x) ~isempty(regexp(x,[str{a},'\_\d+'],'once')), names));
    names_all = [names_all,names_temp];
end
names = names_all;

% Load all optimization results
AllResults = cell(1,length(names));
for a = 1:length(AllResults)
    if ~isempty(names{a})
        strt = load(names{a});
        AllResults{a} = strt.OptimizationResults;
    end
end
AllResults = AllResults(~cellfun(@isempty,AllResults));

% Combine Data
combResults = struct;
for a = 1:length(AllResults)
    fdNames = fieldnames(AllResults{a});
    optstruct = AllResults{a};

    if ~ismember('LacNAcLenPenalty',fdNames)
        LacNAcLenPenalty = [];
    end

    for b = 1:length(fdNames)

        if ~ismember(fdNames{b},fieldnames(combResults))
            combResults.(fdNames{b}) = optstruct.(fdNames{b});
        else
            combResults.(fdNames{b}).xval = [combResults.(fdNames{b}).xval;optstruct.(fdNames{b}).xval];
            combResults.(fdNames{b}).fval = [combResults.(fdNames{b}).fval;optstruct.(fdNames{b}).fval];
            combResults.(fdNames{b}).LacNAcLenPenalty = [combResults.(fdNames{b}).LacNAcLenPenalty;optstruct.(fdNames{b}).LacNAcLenPenalty];

        end
    end
end

end
