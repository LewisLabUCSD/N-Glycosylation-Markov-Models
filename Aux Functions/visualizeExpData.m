function visualizeExpData(DataSet,ProfNames,numSel)

if isempty(ProfNames)
    ProfNames = DataSet.ProfNames;
end


if ischar(ProfNames)
    ProfNames = {ProfNames};
end

profiles = DataSet.profiles;
mz_all = DataSet.mz_all;

if ~isempty(DataSet.LinkageResStruct)
    LinkageResStruct = DataSet.LinkageResStruct;
    LinkageInfoAvail = DataSet.LinkageInfoAvail;
    LinkageResStructSel = DataSet.LinkageResStructSel;
end

for a = 1:length(ProfNames)

    ExpData = profiles(:,strcmp(DataSet.ProfNames,ProfNames{a}));

    if isempty(numSel)
        numSel = sum(ExpData~=0);
    end


    [~,selIdx] = sort(ExpData);
    selIdx = selIdx((end-numSel+1):end);
    [~,ord] = sort(DataSet.mz_all(selIdx));
    selIdx = selIdx(ord);
    ExpData = ExpData(selIdx);
    mz_all = DataSet.mz_all(selIdx);


    figure;
    bar(ExpData,'EdgeColor',[0.25 0.25 0.25],'FaceColor',[0.25 0.25 0.25]);
    xticks(1:length(mz_all));
    xticklabels(round(mz_all));
    title({'Experimental Profile', ['(',strrep(ProfNames{a},'_','/'),')']});
    xtickangle(45);

end

end