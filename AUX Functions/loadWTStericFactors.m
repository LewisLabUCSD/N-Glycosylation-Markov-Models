function WTSteric = loadWTStericFactors(str)


% Identify optimization result files from the Data folder
listing = dir('**/Data/OptimizationResults');
names = {listing.name};
names = names(contains(names,str));

% Load all optimization results
AllResults = cell(1,length(names));
for a = 1:length(AllResults)
    strt = load(names{a});
    AllResults{a} = strt.OptimizationResults;
end

% Combine Data
combResults = struct;
for a = 1:length(AllResults)
   fdNames = fieldnames(AllResults{a});
   optstruct = AllResults{a};
   for b = 1:length(fdNames)  
       if ~ismember(fdNames{b},fieldnames(combResults))
       combResults.(fdNames{b}) = optstruct.(fdNames{b});
       else
           combResults.(fdNames{b}).xval = [combResults.(fdNames{b}).xval;optstruct.(fdNames{b}).xval];
           combResults.(fdNames{b}).fval = [combResults.(fdNames{b}).fval,optstruct.(fdNames{b}).fval];
       end
   end
end

WTSteric = mean(combResults.WT.xval(find(~isoutlier(combResults.WT.fval)),end-length(combResults.WT.OptimizationProblem.stericRxns)+1:end),1);

end