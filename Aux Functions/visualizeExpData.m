function visualizeExpData(DataSet)

ProfNames = DataSet.ProfNames;
profiles = DataSet.profiles;
mz_all = DataSet.mz_all;

if ~isempty(DataSet.LinkageResStruct)
    LinkageResStruct = DataSet.LinkageResStruct;
    LinkageInfoAvail = DataSet.LinkageInfoAvail;
    LinkageResStructSel = DataSet.LinkageResStructSel;
end

for a = 1:length(ProfNames)
    
    ExpData = profiles(:,a);
    figure;
    bar(ExpData);
    xticks(1:length(mz_all));
    xticklabels(round(mz_all));
    title({'Experimental Profile', ['(',ProfNames{a},')']});
    xtickangle(45);

end

end