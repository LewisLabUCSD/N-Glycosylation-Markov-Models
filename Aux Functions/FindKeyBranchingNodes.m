function [AllEdges, AllNodes] = FindKeyBranchingNodes(G, Absnodes)

% Inverse network
Edges = G.Edges.EndNodes(:,[2,1]);
G = digraph(Edges(:,1),Edges(:,2));
NodeSel = G.Nodes.Name;

% Sequentially remove nodes and check if the graph remains connnected
for a = 1:length(NodeSel)
    
    if  ~strcmp(NodeSel{a},'(Ma2Ma2Ma3(Ma2Ma3(Ma2Ma6)Ma6)Mb4GNb4GN);Asn[cg]') && ~ismember(NodeSel{a},Absnodes)
        % gaps to be fused upon removing nodes
        inNodes = Edges(strcmp(Edges(:,2),NodeSel{a}),1);
        outNodes = Edges(strcmp(Edges(:,1),NodeSel{a}),2);
        
        % try removing a node & check connectedness of the graph
        H = graph(G.Edges.EndNodes(:,1),G.Edges.EndNodes(:,2));
        H = rmnode(H,NodeSel{a});
        
        % modify G if appropriate
        if all(conncomp(H) == 1) || (length(inNodes) == 1 && length(outNodes) == 1)
            % remove the node
            G = rmnode(G,NodeSel{a});
            % fuse gaps
            if ~isempty(inNodes) && ~isempty(outNodes)
                for k1 = 1:length(inNodes)
                    for k2 = 1:length(outNodes)
                        G = addedge(G,inNodes(k1),outNodes(k2));
                    end
                end
            end
        end
    end
    
end

% remove no exit glycans
Node_noExt = setdiff(G.Edges.EndNodes(:,2),union(G.Edges.EndNodes(:,1),Absnodes));
G = rmnode(G,Node_noExt);

% Record traversed edges and nodes
AllEdges = [G.Edges.EndNodes(:,2),G.Edges.EndNodes(:,1)];
AllNodes = G.Nodes.Name;

end