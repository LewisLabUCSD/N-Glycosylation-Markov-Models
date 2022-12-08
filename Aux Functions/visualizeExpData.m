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
else
    error('Linkage annotations are incomplete or not provided for m/z values.');
end

for a = 1:length(ProfNames)

    idx = find(strcmp(ProfNames{a},DataSet.ProfNames),1);
    ExpData = profiles(:,strcmp(DataSet.ProfNames,ProfNames{a}));

    if isempty(numSel)
        numSel = sum(ExpData~=0);
    end


    [~,selIdx] = sort(ExpData);
    if numSel>length(selIdx)
        numSel = length(selIdx);
    end
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
    ylabel('Relative Abundance');
    xlabel('m/z (or numbering)');

    % visualize glycan annotations
    if LinkageInfoAvail
        sel = DataSet.LinkageResStructSel(:,idx);
        glys = DataSet.LinkageResStruct(sel);
        mz_temp = DataSet.mz(sel);
        DrawGlycanStructure(glys,['Exp ',ProfNames{a}],mz_temp);
    end

end

end