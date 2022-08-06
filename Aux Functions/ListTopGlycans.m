function T = ListTopGlycans(Prof, OptimizationResults, GenericNetwork, TopSel)

% load data
Predata_raw = OptimizationResults.(Prof).Predata_raw;
Glys_raw = OptimizationResults.(Prof).Glys_raw;
AbsGlyIdx = GenericNetwork.AbsGlyIdx;
Glys = GenericNetwork.Glys;

% sort and select data
[~,idx] = sort(mean(Predata_raw,1),'descend');
Predata_raw = Predata_raw(:,idx(1:TopSel));
Glys_raw = Glys_raw(idx(1:TopSel));
Predata_raw_mean = mean(Predata_raw,1)';
Predata_raw_std = std(Predata_raw,[],1)';


fprintf('***********************************************************************\n');

% print data to command window
mz = nan(size(Glys_raw));
for a = 1:length(Predata_raw)
    mzflag = cellfun(@(x) any(strcmp(Glys_raw{a},Glys(x))),AbsGlyIdx);
    if any(mzflag)
        mz(a) = GenericNetwork.DataSet.mz_all(mzflag);  
    end
end

T = table(mz,Glys_raw,Predata_raw_mean,Predata_raw_std,'VariableNames',{'m/z','Glycoform','Relative Ratio','St. Dev.'});

disp(T);

fprintf('***********************************************************************\n');

end