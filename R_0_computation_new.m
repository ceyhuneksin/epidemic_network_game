%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%% R_0 computation. 
% Runs the dynamics until initial sick individual heals. 
%
%% Initialization
N = 100; trials = 100; c_1_values = 3; c_2_values = 26;
beta_equilibrium_bounded = zeros(N,1); T = 200;
delta = 0.2;
%% storage variables
x_store_bounded_information = zeros(N,T,c_1_values,c_2_values,trials);
total_sick_bounded_information = zeros(c_1_values,c_2_values,trials);
actions_store_bounded = zeros(N,c_1_values,c_2_values,trials);
initial_sick_individual= zeros(c_1_values,c_2_values,trials);
infection_probability_bounded_store=zeros(N,c_1_values,c_2_values,trials);
store_network = zeros(N,N,c_1_values,c_2_values,trials);
number_infected_by_patient_zero = zeros(c_1_values,c_2_values,trials);
%% Start runs

%% Create contact and information networks
seed = rng;
rng('default');
%% generate scale free network for the whole set of runs
preferential_attachment_network;
contact_network = graph_new; % geometric network with rewire probability 0.
%% Start trials for varying c_2 and beta values
c_0 = 1; c_1 = 0.24;
for trial = 1:trials %100 trial for a given rate and empathy constant
    trial
    for count_a = 1:c_1_values
        rate = 0.1+(count_a-1)*0.1; % set infection rate beta {0.1, 0.2, 03} 
        for count_b = 1:c_2_values
            store_network(:,:,count_a, count_b,trial) = contact_network;
            c_2 = 0+(count_b-1)*0.02;% set empathy constant 0:0.02:0.5 
            %% initialize the run setup
            x_bounded = zeros(N,1);
            infection_probability_bounded = zeros(N,1);
            sick_person = randi(N);
            initial_sick_individual(count_a, count_b,trial) = sick_person;
            x_bounded(sick_person) = 1;
            %% Dynamics
            time = 1;
            x_store_bounded_information(:,time,count_a,count_b,trial) = x_bounded;
            
            while x_bounded(sick_person) == 1 % run until sick individual heals          
                beta_equilibrium_bounded = iterated_elimination(contact_network,x_bounded,N,c_0,c_1,c_2);
            %% Compute probabilities of infection      
            %% This looks for patient zero infecting its neighbors
                for agent = 1:N
                    if x_bounded(agent) == 1;
                        infection_probability_bounded(agent) = 1 - delta;
                        x_bounded(agent) = rand < infection_probability_bounded(agent);
                    elseif x_bounded(agent) == 0 && contact_network(agent,sick_person)==0
                        neighbors = find(contact_network(agent,:)==1);
                        staying_healthy_with_each_interaction = 1;
                        %                 infection_probability_bounded(agent,tt) = 1 - (1 - beta_equilibrium_bounded(agent,tt))^(sum(x_bounded(find(contact_network(agent,:)==1),tt)));
                        for neighboring_agent = 1:numel(neighbors)
                            staying_healthy_with_each_interaction = staying_healthy_with_each_interaction*(1- rate*beta_equilibrium_bounded(agent)*beta_equilibrium_bounded(neighbors(neighboring_agent))* x_bounded(neighbors(neighboring_agent)));
                        end
                        infection_probability_bounded(agent) = 1 - staying_healthy_with_each_interaction;
                        x_bounded(agent) = rand < infection_probability_bounded(agent);
                    elseif x_bounded(agent) == 0 && contact_network(agent,sick_person)==1 % if susceptible is a neighbor of patient zero (sick_person)
                        neighbors = find(contact_network(agent,:)==1); % find neighbors of agent
                        for neighboring_agent = 1:numel(neighbors) % among neighbors of agent check who infects it independently
                            if x_bounded(agent) == 0 
                                % if rand< beta x a_i x a_j x state neigh.
                                x_bounded(agent) = rand < rate*beta_equilibrium_bounded(agent)*beta_equilibrium_bounded(neighbors(neighboring_agent))* x_bounded(neighbors(neighboring_agent));
                                if neighbors(neighboring_agent) == sick_person && x_bounded(agent) == 1 % if infected by patient zero
                                    number_infected_by_patient_zero(count_a,count_b,trial) = number_infected_by_patient_zero(count_a,count_b,trial)+1;
                                end
                            end
                        end                        
                    end
                end                
                time = time+1;
                x_store_bounded_information(:,time,count_a,count_b,trial) = x_bounded;               
            end
            %% Store stuff
            total_sick_bounded_information(count_a,count_b,trial) = sum(x_bounded);
            actions_store_bounded(:,count_a,count_b,trial) = beta_equilibrium_bounded;
            infection_probability_bounded_store(:,count_a,count_b,trial) = infection_probability_bounded;
            count_b;
            %         count_b = count_b+1;
        end
        count_a
        %     count_a = count_a+1;
    end
end