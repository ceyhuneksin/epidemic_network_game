%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%% 
% iterated elimination of strictly dominated strategies
function [ respond ] = iterated_elimination(A,x,N,c_0,c_1,c_2)
%
%%

iteration_max = N;
response = zeros(N,1);
agent_list_steps = cell(1,iteration_max);
even_agent_list_steps = cell(1,ceil(iteration_max/2));
odd_agent_list_steps = cell(1,ceil(iteration_max/2));

iteration = 1;
while iteration <= iteration_max
    iteration_list = [];
    union_list = cell2mat(agent_list_steps);
    jj =1;
    for steps = 2:2:iteration
        even_agent_list_steps{jj} = agent_list_steps{steps};
        jj = jj +1;
    end
    jj =1;
    for steps = 1:2:iteration
        odd_agent_list_steps{jj} = agent_list_steps{steps};
        jj = jj +1;
    end
    even_union_list = cell2mat(even_agent_list_steps);
    odd_union_list = cell2mat(odd_agent_list_steps);
    
    if rem(iteration,2) ==1  % if the update is odd
        for ii = 1:N 
            if x(ii)==0 && ismember(ii,union_list) == 0
                count_remaining_sick_neighbors = 0;
                for neighbor = 1:N
                    if A(ii,neighbor) == 1 && ismember(neighbor,even_union_list) == 0 && x(neighbor) == 1
                        count_remaining_sick_neighbors = count_remaining_sick_neighbors+1;
                    end
                end
                if c_0> c_1*count_remaining_sick_neighbors
                    iteration_list = [iteration_list ii];
                    response(ii) = 1;
                end
            elseif x(ii) == 1 && ismember(ii,union_list) == 0
                count_remaining_healthy_neighbors = 0;
                for neighbor = 1:N
                    if A(ii,neighbor) == 1 && ismember(neighbor,even_union_list) == 0 && x(neighbor) == 0
                        count_remaining_healthy_neighbors = count_remaining_healthy_neighbors+1;
                    end
                end
                if c_0 > c_2*count_remaining_healthy_neighbors
                    iteration_list = [iteration_list ii];
                    response(ii) = 1;
                end
            end
        end
    else  % if the update is even
        for ii = 1:N
            if x(ii)==0  && ismember(ii,union_list) == 0
                count_sick_social_neighbors = 0;
                for neighbor = 1:N
                    if A(ii,neighbor) == 1 && ismember(neighbor,odd_union_list) == 1 && x(neighbor) == 1
                        count_sick_social_neighbors = count_sick_social_neighbors+1;
                    end
                end
                if c_0< c_1*count_sick_social_neighbors
                    iteration_list = [iteration_list ii];
                    response(ii) = 0;
                end
            elseif x(ii) == 1  && ismember(ii,union_list) == 0
                count_healthy_social_neighbors = 0;
                for neighbor = 1:N
                    if A(ii,neighbor) == 1 && ismember(neighbor,odd_union_list) == 1 && x(neighbor) == 0
                        count_healthy_social_neighbors = count_healthy_social_neighbors+1;
                    end
                end
                if c_0<c_2*count_healthy_social_neighbors
                    iteration_list = [iteration_list ii];
                    response(ii) = 0;
                end
            end
        end
    end
    agent_list_steps{iteration} = iteration_list;
    if size(iteration_list)==0
        iteration;
        iteration = iteration_max;
    end
    iteration = iteration+1;
end

union_list = cell2mat(agent_list_steps);
not_in_union_list = find(ismember(1:N,union_list) == 0);

% if equilibrium is not unique, let susceptibles socialize when c_2>c_1
if c_1 > c_2
    for ii = 1:numel(not_in_union_list)
        if x(not_in_union_list(ii)) == 0
            response(not_in_union_list(ii)) = 0;
        else
            response(not_in_union_list(ii)) = 1;
        end
    end
else
    for ii = 1:numel(find(ismember(1:N,union_list) == 0))       
        if x(not_in_union_list(ii)) == 0
            response(not_in_union_list(ii)) = 1;
        elseif x(not_in_union_list(ii)) == 1
            response(not_in_union_list(ii)) = 0;
        end
    end
end

respond = response;
%% compute utilities
utility = zeros(N,1);
utility(:) = c_0 * response - c_1 * diag(response)* diag(1-x)* A * diag(x)* response ...
        - c_2 * diag(response) *  diag(x)* A * diag(1-x)* response;

%
welfare = 0;
welfare = sum(x) + sum(response) - 2* x' * response;


%% Check: switching is not profitable
new_utility = zeros(N,1);
for ii = 1:N
    new_actions = response;
    new_actions(ii) = 1- new_actions(ii);
    new_utility(ii) = c_0 * new_actions(ii) - c_1 * new_actions(ii)* (1-x(ii))* A(ii,:) * diag(x)* new_actions ...
        - c_2 * new_actions(ii) *  x(ii) * A(ii,:) * diag(1-x)* new_actions;
end

 sum(new_utility - utility >= 0);