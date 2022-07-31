function [mzRes, linkagePosRes] = PrepGlycanAnnotationConstraints(Prof,DataSet,GenericNetwork)

% Initiate variables
linkagePosRes = {};
mzRes = [];

% Prepare 
Annotation_idx = find(DataSet.LinkageResStructSel(:,strcmp(DataSet.ProfNames,Prof)));
if ~isempty(Annotation_idx)
    [~,mzRes] = ismember(DataSet.mz(Annotation_idx),DataSet.mz_all);
    linkageRes = strcat(DataSet.LinkageResStruct(Annotation_idx),'[ab]');
    linkagePosRes = cellfun(@(x) find(strcmp(x,GenericNetwork.Glys)),linkageRes,'UniformOutput',0);
end

end