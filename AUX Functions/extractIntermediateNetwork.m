function [AllNodes, AllEdges] = extractIntermediateNetwork(G, nodes)
  
   % Reverse network
   Edges = G.Edges.EndNodes(:,[2,1]);
   G = digraph(Edges(:,1),Edges(:,2));

   % Initiate variables
   edges = cell(length(nodes),1);
   
   % record all traversed edges
   for a = 1:length(nodes)
   T = bfsearch(G,nodes{a},'allevents');
   edges{a} = (T.EdgeIndex)';edges{a} = edges{a}(~isnan(edges{a}));
   end
   edges = unique([edges{:}]);
   
   AllEdges = G.Edges.EndNodes(edges,[2,1]);
   AllNodes = union(AllEdges(:,1),AllEdges(:,2));
   
end