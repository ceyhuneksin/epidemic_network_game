%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%% epidemic dynamics with rational response
% Runs dynamics for 50 trials 
% Each trial generates a new scale free network. 
% Each trial contains runs for 11 c_1 and 11 c_2 values
% For a given c_1 and c_2 value runing horizon is 200 steps
%% Initialization
% Number of individuals, Simulation horizon, number of values of c_1, c_2
N = 100; T = 150; trials = 50; c_1_values = 11; c_2_values = 11;
beta_equilibrium_bounded = zeros(N,T); % beta actions
% select healing rate and infection rate beta = 0.1, 0.2, 0.3 to get the figures
delta = 0.2; rate = 0.1; 
%% storage variables
x_store_bounded_information = zeros(N,T,c_1_values,c_2_values,trials); % stores state
actions_store_bounded = zeros(N,T,c_1_values,c_2_values,trials); % stores actions
welfare = zeros(T,c_1_values,c_2_values,trials); % stores welfare
alpha = 0.5; 
aggregate_utility = zeros(T,c_1_values,c_2_values,trials); 
store_eradication_time = zeros(c_1_values,c_2_values,trials);
infection_probability_bounded_store = zeros(N,T,c_1_values,c_2_values,trials); 
store_network = zeros(N,N,c_1_values,c_2_values,trials);
%% Start runs
seed = rng; % select seed so that you get the same result every time. 
rng('default');
%% create network
% TO create a Erdos-Renyi Network
% erdosRenyi_graph(N,p);
% contact_network = G.Adj;
% TO Create a scale free network
% preferential_attachment_network;
% TO Create a geometric/small-world network
% p_r = 0; %rewire probability % set to zero for no rewiring % geometric network
% SmallWorldify
% contact_network = graph_new;
% count_b = 1;
c_0 = 1;
for trial = 1:trials
    trial % print trial
    % generate a new preferential network for each trial
    preferential_attachment_network;
    contact_network = graph_new; % geometric network with rewire probability 0.
for count_a = 1:c_1_values
    % change c_1 (risk aversion constant) values
    c_1 = 0+(count_a-1)*0.1; 
    for count_b = 1:c_2_values
        % change c_2 (empathy constant constant) values
        c_2 = 0+(count_b-1)*0.05;
        % save the network of this trial
        store_network(:,:,count_a, count_b,trial) = contact_network;
        %% initialize state and probability values
        x_bounded = zeros(N,T); % state values
        infection_probability_bounded = zeros(N,T);  % probability of infection
%         x_bounded(randi(N),1) = 1; % select one random individual to be selected initialy
        x_bounded(:,1) = rand(N,1)<1; % start with everyone sick.
        %% Start the disease dynamics
        for tt = 1:T-1
            %%
            if sum(x_bounded(:,tt))==0 % if disease is eradicated skip the rest
                tt % print out the eradication time
                store_eradication_time(count_a,count_b,trial) = tt; % get the eradication time
                x_bounded(:,tt:T) = 0; % set the rest of the disease states to zero.
                beta_equilibrium_bounded(:,tt:T) = 1; % set the rest of actions to 1.
                welfare(tt:T,count_a,count_b,trial) =  alpha*N;
                aggregate_utility(tt:T,count_a,count_b) = N * c_0;
                break
            else            
            %% Agent equilibrium computation
            rand_sample = rand(N,1);
            %% Compute equilibrium using best response dynamics
            beta_equilibrium_bounded(:,tt) = iterated_elimination(contact_network,x_bounded(:,tt),N,c_0,c_1,c_2);
            x = x_bounded(:,tt); 
            % calculate aggregate utility and welfare
            response = beta_equilibrium_bounded(:,tt);
            aggregate_utility(tt,count_a,count_b) = sum(c_0 * response - c_1 * diag(response)* diag(1-x)* contact_network * diag(x)* response ...
                - c_2 * diag(response) *  diag(x)* contact_network * diag(1-x)* response);            
            welfare(tt,count_a,count_b,trial) = (1-alpha)*sum(x_bounded(:,tt)) + alpha*sum(beta_equilibrium_bounded(:,tt))  - x_bounded(:,tt)' * beta_equilibrium_bounded(:,tt);            
            %% Compute probabilities for social distancing model with equilibrium responses
            for agent = 1:N
                if x_bounded(agent,tt) == 1;
                    infection_probability_bounded(agent,tt) = 1 - delta;
                else
                    neighbors = find(contact_network(agent,:)==1);
                    staying_healthy_with_each_interaction = 1;
                    for neighboring_agent = 1:numel(neighbors)
                        staying_healthy_with_each_interaction = staying_healthy_with_each_interaction*(1- rate*beta_equilibrium_bounded(agent,tt)*beta_equilibrium_bounded(neighbors(neighboring_agent),tt)* x_bounded(neighbors(neighboring_agent),tt));
                    end
                    infection_probability_bounded(agent,tt) = 1 - staying_healthy_with_each_interaction;
                end
            end
            % probability that agents get sick or remain infected
            x_bounded(:,tt+1) = rand_sample<infection_probability_bounded(:,tt);
            end
        end
        %% Store stuff
        x_store_bounded_information(:,:,count_a,count_b,trial) = x_bounded;       
        actions_store_bounded(:,:,count_a,count_b,trial) = beta_equilibrium_bounded;
        infection_probability_bounded_store(:,:,count_a,count_b,trial) = infection_probability_bounded;        
        count_b;
    end
    count_a % print 
end
end