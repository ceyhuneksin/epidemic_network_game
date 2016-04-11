%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%%%%%%% Preferential attachment %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr_nodes =100;
%% generate network
graph_new = spalloc(nr_nodes, nr_nodes, 50);
degree_nodes = zeros(nr_nodes,1);

graph_new(1,2) = 1; % set two initial nodes to be connected to one another
graph_new(2,1) = 1; % set two initial nodes to be connected to one another
degree_nodes(1) = 1;
degree_nodes(2) = 1;
cumulative_degree_distribution = zeros(nr_nodes,1);

cumulative_degree_distribution(1) = degree_nodes(1)/sum(degree_nodes);
for node = 2:nr_nodes
    cumulative_degree_distribution(node) = cumulative_degree_distribution(node-1) + degree_nodes(node)/sum(degree_nodes);
end


% add nodes and connect one at a time 
for node = 3:nr_nodes
    indices = find(cumulative_degree_distribution>rand);
    graph_new(indices(1),node) = 1;
    graph_new(node,indices(1)) = 1;
    degree_nodes(indices(1)) = degree_nodes(indices(1))+1;
    degree_nodes(node) = 1;
    %% update cumulative distribution
    cumulative_degree_distribution(1) = degree_nodes(1)/sum(degree_nodes);
    for node = 2:nr_nodes
        cumulative_degree_distribution(node) = cumulative_degree_distribution(node-1) + degree_nodes(node)/sum(degree_nodes);
    end
end

