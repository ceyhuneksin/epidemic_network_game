%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%% epidemic dynamics game theory
%
%% Initialization
N = 4; T = 20; trials = 1; c_1_values = 1; c_2_values = 2;
beta_equilibrium_bounded = zeros(N,T); beta_worst_case = zeros(N,T); beta_welfare = zeros(N,T); beta_noempathy = zeros(N,T);
beta_best_response = zeros(N,T);
delta = 0.2; rate = 0.4;
%% storage variables
x_store_bounded_information = zeros(N,T,c_1_values,c_2_values,trials); x_bench_store = zeros(N,T,c_1_values,c_2_values,trials); x_store_worst = zeros(N,T,c_1_values,c_2_values,trials); x_store_welfare = zeros(N,T,c_1_values,c_2_values,trials);x_store_noempathy = zeros(N,T,c_1_values,c_2_values,trials);
x_store_best_response = zeros(N,T,c_1_values,c_2_values,trials);
actions_store_bounded = zeros(N,T,c_1_values,c_2_values,trials); actions_store_worst = zeros(N,T,c_1_values,c_2_values,trials);actions_store_welfare = zeros(N,T,c_1_values,c_2_values,trials);actions_store_noempathy = zeros(N,T,c_1_values,c_2_values,trials);
actions_store_best_response = zeros(N,T,c_1_values,c_2_values,trials);
welfare = zeros(T,c_1_values,c_2_values,trials); welfare_bench = zeros(T,c_1_values,c_2_values,trials); welfare_worst = zeros(T,c_1_values,c_2_values,trials); welfare_welfare = zeros(T,c_1_values,c_2_values,trials); welfare_noempathy = zeros(T,c_1_values,c_2_values,trials);
welfare_best_response = zeros(T,c_1_values,c_2_values,trials);
alpha = 0.5;
aggregate_utility = zeros(T,c_1_values,c_2_values,trials); aggregate_utility_bench = zeros(T,c_1_values,c_2_values,trials);aggregate_utility_worst = zeros(T,c_1_values,c_2_values,trials); aggregate_utility_welfare = zeros(T,c_1_values,c_2_values,trials); aggregate_utility_noempathy = zeros(T,c_1_values,c_2_values,trials);
aggregate_utility_best_response = zeros(T,c_1_values,c_2_values,trials);
store_eradication_time = zeros(c_1_values,c_2_values,trials);

infection_probability_bounded_store=zeros(N,T,c_1_values,c_2_values,trials); infection_probability_bench_store =zeros(N,T,c_1_values,c_2_values,trials); infection_probability_worst_store=zeros(N,T,c_1_values,c_2_values,trials); infection_probability_welfare_store=zeros(N,T,c_1_values,c_2_values,trials);infection_probability_noempathy_store=zeros(N,T,c_1_values,c_2_values,trials);
infection_probability_best_response_store=zeros(N,T,c_1_values,c_2_values,trials);
store_network = zeros(N,N,c_1_values,c_2_values,trials);
%% Create contact networks for the demo

contact_network = [0 1 1 1; 1 0 0 0; 1 0 0 0; 1 0 0 0]; % star network
%% Start runs
rand('state',15)
c_0 = 1;
for trial = 1:trials
    trial
for count_a = 1:c_1_values
    c_1 = 0.4;
%     rate = count_a/10;
%     c_1 = 0.21;
    for count_b = 1:c_2_values
        store_network(:,:,count_a, count_b,trial) = contact_network;
        c_2 = 0.2+(count_b-1)*0.2;
        %%
        x_bounded = zeros(N,T); 
        infection_probability_bounded = zeros(N,T);
        x_bounded(1,1) = 1; % start with center node sick.
        x_bench(:,1) = x_bounded(:,1);
        %% Dynamics
        for tt = 1:T-1
            %%
            if sum(x_bounded(:,tt))==0 % disease is eradicated skip the rest
                tt
                store_eradication_time(count_a,count_b,trial) = tt;
                x_bounded(:,tt:T) = 0;
                beta_equilibrium_bounded(:,tt:T) = 1;
                welfare(tt:T,count_a,count_b,trial) =  alpha*N;
                aggregate_utility(tt:T,count_a,count_b) = N * c_0;
                break
            else            
            %% Agent equilibrium computation
            rand_sample = rand(N,1);
            %% Compute equilibrium using best response dynamics
            beta_equilibrium_bounded(:,tt) = iterated_elimination(contact_network,x_bounded(:,tt),N,c_0,c_1,c_2);

            x = x_bounded(:,tt); response = beta_equilibrium_bounded(:,tt);
            aggregate_utility(tt,count_a,count_b) = sum(c_0 * response - c_1 * diag(response)* diag(1-x)* contact_network * diag(x)* response ...
                - c_2 * diag(response) *  diag(x)* contact_network * diag(1-x)* response);
            welfare(tt,count_a,count_b,trial) = (1-alpha)*sum(x_bounded(:,tt)) + alpha*sum(beta_equilibrium_bounded(:,tt))  - x_bounded(:,tt)' * beta_equilibrium_bounded(:,tt);
            
            %% Compute probabilities for social distancing model with bounded social distancing
            for agent = 1:N
                if x_bounded(agent,tt) == 1;
                    infection_probability_bounded(agent,tt) = 1 - delta;
                else
                    neighbors = find(contact_network(agent,:)==1);
                    staying_healthy_with_each_interaction = 1;
                    %                 infection_probability_bounded(agent,tt) = 1 - (1 - beta_equilibrium_bounded(agent,tt))^(sum(x_bounded(find(contact_network(agent,:)==1),tt)));
                    for neighboring_agent = 1:numel(neighbors)
                        staying_healthy_with_each_interaction = staying_healthy_with_each_interaction*(1- rate*beta_equilibrium_bounded(agent,tt)*beta_equilibrium_bounded(neighbors(neighboring_agent),tt)* x_bounded(neighbors(neighboring_agent),tt));
                    end
                    infection_probability_bounded(agent,tt) = 1 - staying_healthy_with_each_interaction;
                end
            end
            x_bounded(:,tt+1) = rand_sample<infection_probability_bounded(:,tt);

            end
        end
        %% Store stuff
        x_store_bounded_information(:,:,count_a,count_b,trial) = x_bounded;        
        actions_store_bounded(:,:,count_a,count_b,trial) = beta_equilibrium_bounded;
        infection_probability_bounded_store(:,:,count_a,count_b,trial) = infection_probability_bounded;
        count_b;
    end
    count_a
end
end