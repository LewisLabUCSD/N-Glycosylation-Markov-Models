function path = ConcentrationGradientChaser(G, NodeWt, GlySel)

% reverse network
Edges_old = G.Edges.EndNodes(:,[2,1]);
Nodes_old = G.Nodes.Name;
G = digraph(Edges_old(:,1),Edges_old(:,2));
Edges = G.Edges.EndNodes;
Nodes = G.Nodes.Name;
NodeWt = cellfun(@(x) NodeWt(strcmp(x,Nodes_old)), Nodes);

% Initiate variable for concentration gradient chasing
CurrNode = GlySel;
continueFlag = true;
path = {GlySel};

% start chasing
while continueFlag
    
    NextNode_sel = find(strcmp(CurrNode,Edges(:,1)));
    
    if isempty(NextNode_sel)
        continueFlag = false;
    else
    NextNode = Edges(NextNode_sel,2);
    NextNode_wt = cellfun(@(x) NodeWt(strcmp(x,Nodes)), NextNode);
    
    % chase the path crossing the nodes with the biggest
    % pseudo-concentrations
    [~,idx_max] = max(NextNode_wt);
    
    % record & update
    path = [path,NextNode(idx_max)];
    CurrNode = NextNode{idx_max};
    end    
    
end

% reverse path
path = path(end:-1:1);

end