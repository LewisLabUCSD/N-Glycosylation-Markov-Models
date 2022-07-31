function MajorEdges = RemoveNonBranchingNodes(MajorEdges)

nodes = unique([MajorEdges(:,1);MajorEdges(:,2)]);
nonBranchingNodes = cellfun(@(x) sum(strcmp(x,MajorEdges(:,1))) == 1 && sum(strcmp(x,MajorEdges(:,2))) == 1,nodes );

% remove one non-branching node each time and fuse the gap
while ~all(nonBranchingNodes==0)
    
    idx = find(nonBranchingNodes,1);
    
    % inNode
    inNodeRow = strcmp(MajorEdges(:,2),nodes(idx));
    inNode = MajorEdges{inNodeRow,1};

    % outNode
    outNodeRow = strcmp(MajorEdges(:,1),nodes(idx));
    outNode = MajorEdges{outNodeRow,2};
    
    % remove connecting edged
    MajorEdges([find(inNodeRow,1),find(outNodeRow,1)],:) = [];
    
    % fuse gap
    MajorEdges = [MajorEdges;{inNode,outNode}];
    
    % recheck branching point 
    nodes = unique([MajorEdges(:,1),MajorEdges(:,2)]);
    nonBranchingNodes = cellfun(@(x) sum(strcmp(x,MajorEdges(:,1))) == 1 && sum(strcmp(x,MajorEdges(:,2))) == 1,nodes );
end
end
